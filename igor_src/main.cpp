/*
 * main.cpp
 *
 *  Created on: Jan 12, 2015
 *      Author: Quentin Marcou
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "../config.h"
#include "Deletion.h"
#include "Insertion.h"
#include "Genechoice.h"
#include "Model_Parms.h"
#include "Rec_Event.h"
#include "Singleerrorrate.h"
#include "Model_marginals.h"
#include <iostream>
#include "Aligner.h"
#include "GenModel.h"
#include "Dinuclmarkov.h"
#include "Counter.h"
#include "Coverageerrcounter.h"
#include "Bestscenarioscounter.h"
#include "Pgencounter.h"
#include "Errorscounter.h"
#include "Utils.h"
#include <chrono>
#include <set>

#include <string>
#include "CDR3SeqData.h"
#include "ExtractFeatures.h"

using namespace std;


// FIXME: PUT THIS IN ANOTHER FILE


// TODO: Possible typedef definitions for code readability.
typedef std::string strSeqID;    // fasta description >strSeqID
typedef std::string strSequence;
typedef vector< pair<const int, const strSequence> > VectorIndexedSeq;
typedef pair<strSeqID, strSequence> PairSeq; // FIXME: Not a good name convention.
typedef vector<PairSeq>  VectorGenomicTemplate;


int terminate_IGoR_with_error_message(const forward_list<string> error_messages){
	string igor_error_prefix = "[IGoR] ERROR: ";
	for(string error_message : error_messages){
		cerr<<igor_error_prefix<<error_message<<endl;
	}
	cerr<<igor_error_prefix<<"Use \"man igor\", \"igor -help\" or visit "<<PACKAGE_URL<<" to see available commands and their effects."<<endl;
	cerr<<igor_error_prefix<<"Please report any bug by opening an issue on "<<PACKAGE_URL<<" or email: "<<PACKAGE_BUGREPORT<<endl;
	cerr<<igor_error_prefix<<"Terminating IGoR..."<<endl;
	return EXIT_FAILURE;
}

int terminate_IGoR_with_error_message(string error_message){
	return terminate_IGoR_with_error_message(forward_list<string>(1,error_message));
}

// Output both the exception handling message and the actual exception message
int terminate_IGoR_with_error_message(string error_message, exception& e){
	forward_list<string> error_messages;
	error_messages.emplace_front(e.what());
	error_messages.emplace_front(error_message);
	return terminate_IGoR_with_error_message(error_messages);
}

// Output both the exception handling message and the actual exception message
int terminate_IGoR_with_error_message(forward_list<string> error_messages, exception& e){
	error_messages.emplace_front(e.what());
	return terminate_IGoR_with_error_message(error_messages);
}

int main(int argc , char* argv[]){

	//Command line argument iterator
	size_t carg_i = 1;
	//cout<<argv[argc-1]<<endl;

	// Test if some commands were supplied to IGoR
	if(argc<2){
		return terminate_IGoR_with_error_message("The user did not supply IGoR any command.");
	}

	//Task variables
	bool run_demo = false;
	bool align = false;
	bool infer = false;
	bool evaluate = false;
	bool generate = false;
	bool custom = false;

	//Common vars
	string batchname="";


	//Working directory vars
	bool wd = false;
	string cl_path;

	//Chains vars
	bool chain_provided = false;
	string chain_arg_str;
	string chain_path_str;
	bool has_D = false;

	//Species vars
	bool species_provided = false;
	string species_str = "";

	//Custom genomic templates loading variables
	bool custom_v = false;
	string custom_v_path;
	bool custom_d = false;
	string custom_d_path;
	bool custom_j = false;
	string custom_j_path;

	//Custom gene anchors loading variables
	bool custom_v_anchors = false;
	string custom_v_anchors_path;
	bool custom_j_anchors = false;
	string custom_j_anchors_path;

	//Input sequences variables
	bool read_seqs = false;
	Fileformat seqs_fileformat;
	string input_seqs_file;

	//Genomic templates list and aligns parms
	vector<pair<string,string>> v_genomic;
	vector<pair<string,string>> d_genomic;
	vector<pair<string,string>> j_genomic;
	unordered_map<string,size_t> v_CDR3_anchors;
	unordered_map<string,size_t> j_CDR3_anchors;
	bool align_data_is_CDR3 = false;

	//Model parms and marginals
	bool load_last_inferred_parms = false;
	bool custom_cl_parms = false;
	Model_Parms cl_model_parms;
	Model_marginals cl_model_marginals;
	map<size_t,shared_ptr<Counter>> cl_counters_list;

	//Sequence generation parms
	size_t generate_n_seq;
	bool generate_werr = true;
	bool gen_output_CDR3_data = false;
	string gen_filename_prefix="";
	int gen_random_engine_seed=-1; //-1 will cause IGoR to generate a time stamp based seed

	//Inference parms
	bool viterbi_inference = false;
	double likelihood_thresh_inference = 1e-60;
	double proba_threshold_ratio_inference = 1e-5;
	size_t n_iter_inference = 5;
	bool infer_only = false;
	bool no_infer = false;
	set<string> infer_restrict_nicknames;
	bool fix_err_rate = false;
	bool subsample_seqs = false;
	size_t n_subsample_seqs;

	//Sequence evaluation parms
	bool viterbi_evaluate = false;
	double likelihood_thresh_evaluate = 1e-60;;
	double proba_threshold_ratio_evaluate = 1e-5;

	//Alignment parameters
	double heavy_pen_nuc44_vect [] = { // A,C,G,T,R,Y,K,M,S,W,B,D,H,V,N
	        5,-14,-14,-14,-14,2,-14,2,2,-14,-14,1,1,1,0,
	        -14,5,-14,-14,-14,2,2,-14,-14,2,1,-14,1,1,0,
	        -14,-14,5,-14,2,-14,2,-14,2,-14,1,1,-14,1,0,
	        -14,-14,-14,5,2,-14,-14,2,-14,2,1,1,1,-14,0,
	        -14,-14,2,2,1.5,-14,-12,-12,-12,-12,1,1,-13,-13,0,
	        2,2,-14,-14,-14,1.5,-12,-12,-12,-12,-13,-13,1,1,0,
	        -14,2,2,-14,-12,-12,1.5,-14,-12,-12,1,-13,-13,1,0,
	        2,-14,-14,2,-12,-12,-14,1.5,-12,-12,-13,1,1,-13,0,
	        2,-14,2,-14,-12,-12,-12,-12,1.5,-14,-13,1,-13,1,0,
	        -14,2,-14,2,-12,-12,-12,-12,-14,1.5,1,-13,1,-13,0,
	        -14,1,1,1,1,-13,1,-13,-13,1,0.5,-12,-12,-12,0,
	        1,-14,1,1,1,-13,-13,1,1,-13,-12,0.5,-12,-12,0,
	        1,1,-14,1,-13,1,-13,1,-13,1,-12,-12,0.5,-12,0,
	        1,1,1,-14,-13,1,1,-13,1,-13,-12,-12,-12,0.5,0,
	        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	Matrix<double> heavy_pen_nuc44_sub_matrix(15,15,heavy_pen_nuc44_vect);

		//V alignment vars
		bool align_v = false;
		string v_align_filename = "V_alignments.csv";
		double v_align_thresh_value = 50.0;
		Matrix<double> v_subst_matrix = heavy_pen_nuc44_sub_matrix;
		double v_gap_penalty = 50.0;
		bool v_best_align_only = true;
		bool v_best_gene_only = false;
		int v_left_offset_bound = INT16_MIN;
		int v_right_offset_bound = INT16_MAX;
		unordered_map<string,pair<int,int>> v_template_bounds_map;
		bool v_reversed_offsets = false;

		//D alignment vars
		bool align_d = false;
		string d_align_filename = "D_alignments.csv";
		double d_align_thresh_value = 15.0;
		Matrix<double> d_subst_matrix = heavy_pen_nuc44_sub_matrix;
		double d_gap_penalty = 50.0;
		bool d_best_align_only = false;
		bool d_best_gene_only = false;
		int d_left_offset_bound = INT16_MIN;
		int d_right_offset_bound = INT16_MAX;
		unordered_map<string,pair<int,int>> d_template_bounds_map;
		bool d_reversed_offsets = false;

		//J alignment vars
		bool align_j = false;
		string j_align_filename = "J_alignments.csv";
		double j_align_thresh_value = 15.0;
		Matrix<double> j_subst_matrix = heavy_pen_nuc44_sub_matrix;
		double j_gap_penalty = 50.0;
		bool j_best_align_only = true;
		bool j_best_gene_only = false;
		int j_left_offset_bound = INT16_MIN;
		int j_right_offset_bound = INT16_MAX;
		unordered_map<string,pair<int,int>> j_template_bounds_map;
		bool j_reversed_offsets = false;
		
		// Flag to extract CDR3 from aligned sequences.
		bool bFeature_CDR3 = true;


	while(carg_i<argc){

		//Command line argument asking for help
		if(string(argv[carg_i]) == string("-h")
				or string(argv[carg_i]) == string("-help")){
			// Show the manual using man and installed manpages
			system(&("man igor")[0]);
			// End the program without error
			return 0;
		}

		//Command line argument asking for help
		if(string(argv[carg_i]) == string("-v")
				or string(argv[carg_i]) == string("-version")){
			// Display IGoR's version
			cout<<"IGoR version "<<PACKAGE_VERSION<<endl;
			clog<<"Visit "<<PACKAGE_URL<<" to check the latest version!"<<endl;
			// End the program without error
			return 0;
		}

		//Command line argument setting the number of threads
		if(string(argv[carg_i]) == string("-threads")){
			omp_set_num_threads(std::stoi(argv[++carg_i]));
			clog<<"Setting number of threads to: "<<argv[carg_i]<<endl;
		}

		//Command line to redirect the standard output to a file
		else if(string(argv[carg_i]) == string("-stdout_f")){
			clog<<"Redirecting output to file: "<<argv[carg_i+1]<<endl;
			freopen(argv[++carg_i],"a+",stdout);
		}

		else if(string(argv[carg_i]) == string("-batch")){
			++carg_i;
			batchname = argv[carg_i];
			if(batchname[batchname.size()-1] != '_'){
				batchname.append("_");
			}
			clog<<"Batch name set to: "<<batchname<<endl;
		}

		//Set custom genomic template
		else if(string(argv[carg_i]) == string("-set_genomic")){
			//Throw an error if not found?
			while((carg_i+1<argc)
					and (string(argv[carg_i+1]).size()>2)
					and string(argv[carg_i+1]).substr(0,2) == "--"){

				++carg_i;

				if(string(argv[carg_i]) == "--V"){
					custom_v = true;
					++carg_i;
					custom_v_path = string(argv[carg_i]);
				}
				else if(string(argv[carg_i]) == "--D"){
					custom_d = true;
					++carg_i;
					custom_d_path = string(argv[carg_i]);
				}
				else if(string(argv[carg_i]) == "--J"){
					custom_j = true;
					++carg_i;
					custom_j_path = string(argv[carg_i]);
				}
				else{
					return terminate_IGoR_with_error_message("Unknown gene argument \"" + string(argv[carg_i]) +"\" to set genomic templates");
				}
			}

			if((not custom_v)
					and (not custom_d)
					and (not custom_j)){
				return terminate_IGoR_with_error_message("No gene argument was passed after -set_genomic");
			}
		}

		else if(string(argv[carg_i]) == string("-set_CDR3_anchors")){
			while((carg_i+1<argc)
					and (string(argv[carg_i+1]).size()>2)
					and string(argv[carg_i+1]).substr(0,2) == "--"){

				++carg_i;

				if(string(argv[carg_i]) == "--V"){
					custom_v_anchors = true;
					++carg_i;
					custom_v_anchors_path = string(argv[carg_i]);
				}
				else if(string(argv[carg_i]) == "--J"){
					custom_j_anchors = true;
					++carg_i;
					custom_j_anchors_path = string(argv[carg_i]);
				}
				else{
					return terminate_IGoR_with_error_message("Unknown gene argument \"" + string(argv[carg_i]) +"\" to set CDR3 anchors");
				}
			}

			if((not custom_v_anchors)
					and (not custom_j_anchors)){
				return terminate_IGoR_with_error_message("No gene argument was passed after -set_CDR3_anchors");
			}
		}

		else if(string(argv[carg_i]) == "-set_custom_model"){
			custom_cl_parms = true;
			++carg_i;
			try{
				cl_model_parms.read_model_parms(string(argv[carg_i]));
			}catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading custom model parms after \"-set_custom_model\"",  e);
			}
			cl_model_marginals = Model_marginals(cl_model_parms);
			if((carg_i+1)<argc
					and (string(argv[carg_i+1]).substr(0,1)!=string("-"))){
				++carg_i;
				//Check if the next argument is a new command or the corresponding marginal file
				try{
					cl_model_marginals.txt2marginals(string(argv[carg_i]),cl_model_parms);
				}
				catch(exception& e){
					forward_list<string> error_messages;
					error_messages.emplace_front("If you have altered the default model structure corresponding marginals will be created if no model marginals file is passed.");
					error_messages.emplace_front("Make sure file exist or that supplied Model_Parms and Model_Marginals match.");
					error_messages.emplace_front("Exception caught while reading custom marginals after \"-set_custom_model\"");
					return terminate_IGoR_with_error_message(error_messages, e);
				}
			}
			else{
				clog<<"No model marginals file was provided with the custom model parameters, initializing corresponding marginals to a uniform distribution!"<<endl;
				cl_model_marginals.uniform_initialize(cl_model_parms);
			}

			//Check if the model contains a D gene event in order to load the alignments
			auto events_map = cl_model_parms.get_events_map();
			if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side))>0){
				has_D = true;
			}
		}

		else if(string(argv[carg_i]) == "-load_last_inferred"){
			//Cannot read straight here since the working directory have not been defined yet
			load_last_inferred_parms = true;
		}

		else if(string(argv[carg_i]) == "-run_demo"){
			clog<<"Running demo code"<<endl;
			run_demo = true;
		}

		else if(string(argv[carg_i]) == "-run_custom"){
			clog<<"Running custom code"<<endl;
			custom = true;
		}

		else if(string(argv[carg_i]) == "-subsample"){
			subsample_seqs = true;
			++carg_i;
			try{
				n_subsample_seqs = stoi(argv[carg_i]);
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Expected an integer for the number of sequences to subsample, received: \"" + string(argv[carg_i]) + "\"");
			}
		}

		else if(string(argv[carg_i]) == "-align"){
			//Provide a boolean for aligning
			align = true;
			//Check for additional parameters specific to each gene
			while(	(carg_i+1<argc)
					and (string(argv[carg_i+1]).size()>2)
					and (string(argv[carg_i+1]).substr(0,2) == "--")){
				++carg_i;
				string gene_str_val = string(argv[carg_i]);

				//Define init booleans and corresponding variables
				bool thresh_provided = false;
				double thresh_value;

				bool matrix_provided = false;
				//Matrix subst_matrix;

				bool gap_penalty_provided = false;
				double gap_penalty;

				bool best_align_only_provided = false;
				bool best_gene_only_provided = false;
				bool best_align_only;
				bool best_gene_only;

				bool offset_bounds_provided = false;
				bool template_bounds_provided = false;
				int left_offset_bound;
				int right_offset_bound;
				bool reversed_offset_provided = false;
				bool reversed_offsets;
				unordered_map<string,pair<int,int>> template_bounds_map;

				if( (gene_str_val == "--V")
						or (gene_str_val == "--D")
						or (gene_str_val == "--J")
						or (gene_str_val == "--all")){

					while( (carg_i+1<argc)
							and (string(argv[carg_i+1]).size()>3)
							and (string(argv[carg_i+1]).substr(0,3) == string("---"))){

						++carg_i;

						if(string(argv[carg_i]) == "---thresh"){
							//Read the alignment score threshold
							++carg_i;
							try{
								thresh_value = stod(string(argv[carg_i]));
							}
							catch (exception& e) {
								return terminate_IGoR_with_error_message("Expected a float for the alignment score threshold, received: \"" + string(argv[carg_i]) + "\"");
							}
							thresh_provided = true;
						}
						else if(string(argv[carg_i]) == "---matrix"){
							//Read the substitution matrix from a file
							++carg_i;
							string matrix_filename = string(argv[carg_i]);
							try{
								if((gene_str_val == "--V") or (gene_str_val == "--all")){
									v_subst_matrix = read_substitution_matrix(matrix_filename);
								}
								if((gene_str_val == "--D") or (gene_str_val == "--all")){
									d_subst_matrix = read_substitution_matrix(matrix_filename);
								}
								if((gene_str_val == "--J") or (gene_str_val == "--all")){
									j_subst_matrix = read_substitution_matrix(matrix_filename);
								}
							}
							catch(exception& e){
								return terminate_IGoR_with_error_message("Exception caught while reading the provided substitution matrix for gene: " + gene_str_val,  e);
							}

						}
						else if(string(argv[carg_i]) == "---gap_penalty"){
							//Read the gap penalty
							++carg_i;
							try{
								gap_penalty = stod(string(argv[carg_i]));
							}
							catch (exception& e) {
								return terminate_IGoR_with_error_message("Expected a float for the alignment gap penalty, received: \"" + string(argv[carg_i]) + "\"");
							}
							gap_penalty_provided = true;
						}
						else if(string(argv[carg_i]) == "---best_align_only"){
							//Set the best alignment only boolean
							++carg_i;
							if(string(argv[carg_i]) == "true"){
								best_align_only = true;
							}
							else if(string(argv[carg_i]) == "false"){
								best_align_only = false;
							}
							else{
								return terminate_IGoR_with_error_message("Unkown argument received\"" + string(argv[carg_i]) + "\"to set best alignment only boolean, existing values are: true, false");
							}
							best_align_only_provided = true;

						}
						else if(string(argv[carg_i]) == "---best_gene_only"){
							//Set the best alignment only boolean
							++carg_i;
							if(string(argv[carg_i]) == "true"){
								best_gene_only = true;
							}
							else if(string(argv[carg_i]) == "false"){
								best_gene_only = false;
							}
							else{
								return terminate_IGoR_with_error_message("Unkown argument received\"" + string(argv[carg_i]) + "\"to set best gene/allele candidate only boolean, existing values are: true, false");
							}
							best_gene_only_provided = true;

						}
						else if(string(argv[carg_i]) == "---template_spec_offset_bounds"){
							//Read the offset bounds
							++carg_i;
							try{
								template_bounds_map = read_template_specific_offset_csv(string(argv[carg_i]));
							}
							catch (exception& e) {
								return terminate_IGoR_with_error_message("Exception caught reading a semi-colon separated offset file \"" + string(argv[carg_i]) + "\" for template specific offset bounds.",e);
							}

							template_bounds_provided = true;
						}
						else if(string(argv[carg_i]) == "---offset_bounds"){
							//Read the offset bounds
							++carg_i;
							try{
								left_offset_bound = stoi(string(argv[carg_i]));
							}
							catch (exception& e) {
								return terminate_IGoR_with_error_message("Expected an integer for the left offset bound, received: \"" + string(argv[carg_i]) + "\"");
							}

							++carg_i;
							try{
								right_offset_bound = stoi(string(argv[carg_i]));
							}
							catch (exception& e) {
								return terminate_IGoR_with_error_message("Expected an integer for the right offset bound, received: \"" + string(argv[carg_i]) + "\"");
							}
							offset_bounds_provided = true;
						}
						else if(string(argv[carg_i]) == "---reversed_offsets"){
							//Set the best alignment only boolean
							++carg_i;
							if(string(argv[carg_i]) == "true"){
								reversed_offsets = true;
							}
							else if(string(argv[carg_i]) == "false"){
								reversed_offsets = false;
							}
							else{
								return terminate_IGoR_with_error_message("Unkown argument received\"" + string(argv[carg_i]) + "\"to set reversed offsets boolean, existing values are: true, false");
							}
							reversed_offset_provided = true;
						}
						else{
							return terminate_IGoR_with_error_message("Unknown parameter\"" + string(argv[carg_i]) +"\" for gene " + gene_str_val + " in -align " );
						}
					}
				}
				else if(string(argv[carg_i]) == "--ntCDR3"){
					align_data_is_CDR3 = true;
					// Set the the V and J align thresholds accordingly
					v_align_thresh_value = 0;
					j_align_thresh_value = 0;
				}
				else{
					return terminate_IGoR_with_error_message("Unknown gene specification\"" + string(argv[carg_i]) + "\"for -align");
				}

				//Now assign back the values to the correct variables
				if( (gene_str_val == "--V") or (gene_str_val == "--all")){
					align_v = true;
					if(thresh_provided){
						v_align_thresh_value = thresh_value;
					}

					if(matrix_provided){
						//v_subst_matrix = subst_matrix;
					}

					if(gap_penalty_provided){
						v_gap_penalty = gap_penalty;
					}

					if(best_align_only_provided){
						v_best_align_only = best_align_only;
					}

					if(best_gene_only_provided){
						v_best_gene_only = best_gene_only;
					}

					if(template_bounds_provided){
						v_template_bounds_map = template_bounds_map;
					}

					if(offset_bounds_provided){
						v_left_offset_bound = left_offset_bound;
						v_right_offset_bound = right_offset_bound;
					}

					if(reversed_offset_provided){
						v_reversed_offsets = reversed_offsets;
					}
				}
				if( (gene_str_val == "--D") or (gene_str_val == "--all")){
					align_d = true;
					if(thresh_provided){
						d_align_thresh_value = thresh_value;
					}

					if(matrix_provided){
						//d_subst_matrix = subst_matrix;
					}

					if(gap_penalty_provided){
						d_gap_penalty = gap_penalty;
					}

					if(best_align_only_provided){
						d_best_align_only = best_align_only;
					}

					if(best_gene_only_provided){
						d_best_gene_only = best_gene_only;
					}

					if(template_bounds_provided){
						d_template_bounds_map = template_bounds_map;
					}

					if(offset_bounds_provided){
						d_left_offset_bound = left_offset_bound;
						d_right_offset_bound = right_offset_bound;
					}

					if(reversed_offset_provided){
						d_reversed_offsets = reversed_offsets;
					}
				}
				if( (gene_str_val == "--J") or (gene_str_val == "--all")){
					align_j = true;
					if(thresh_provided){
						j_align_thresh_value = thresh_value;
					}

					if(matrix_provided){
						//j_subst_matrix = subst_matrix;
					}

					if(gap_penalty_provided){
						j_gap_penalty = gap_penalty;
					}

					if(best_align_only_provided){
						j_best_align_only = best_align_only;
					}

					if(best_gene_only_provided){
						j_best_gene_only = best_gene_only;
					}

					if(template_bounds_provided){
						j_template_bounds_map = template_bounds_map;
					}

					if(offset_bounds_provided){
						j_left_offset_bound = left_offset_bound;
						j_right_offset_bound = right_offset_bound;
					}

					if(reversed_offset_provided){
						j_reversed_offsets = reversed_offsets;
					}
				}

			}
		}
 
		else if( (string(argv[carg_i]) == "-infer") or (string(argv[carg_i]) == "-evaluate")){
			//Provide a boolean for inference
			if(string(argv[carg_i]) == "-infer"){
				infer = true;
			}
			else{
				evaluate = true;
			}

			while( (carg_i+1<argc)
					and (string(argv[carg_i+1]).size()>2)
					and string(argv[carg_i+1]).substr(0,2)=="--"){
				++carg_i;
				//Some inference parameters are passed
				if(string(argv[carg_i]) == "--L_thresh"){
					double l_thresh;
					++carg_i;
					try{
						l_thresh = stod(string(argv[carg_i]));
					}
					catch(exception& e){
						return terminate_IGoR_with_error_message("Expected a float for the likelihood threshold, received: \"" + string(argv[carg_i]) + "\"");
					}

					if(infer){
						likelihood_thresh_inference = l_thresh;
					}
					else{
						likelihood_thresh_evaluate = l_thresh;
					}
				}
				else if(string(argv[carg_i]) == "--MLSO"){
					if(infer){
						viterbi_inference = true;
					}
					else{
						viterbi_evaluate = true;
					}
				}
				else if(string(argv[carg_i]) == "--P_ratio_thresh"){
					double p_ratio;
					++carg_i;
					try{
						p_ratio = stod(string(argv[carg_i]));
					}
					catch(exception& e){
						return terminate_IGoR_with_error_message("Expected a float for the probability ratio threshold, received: \"" + string(argv[carg_i]) + "\"");
					}

					if(infer){
						proba_threshold_ratio_inference = p_ratio;
					}
					else{
						proba_threshold_ratio_evaluate = p_ratio;
					}

				}
				else if(string(argv[carg_i]) == "--N_iter"){
					if(infer){
						++carg_i;
						try{
							n_iter_inference = stoi(string(argv[carg_i]));
						}
						catch(exception& e){
							return terminate_IGoR_with_error_message("Expected an integer for the number of iterations to perform for the inference, received: \"" + string(argv[carg_i]) + "\"");
						}
					}
					else{
						return terminate_IGoR_with_error_message("Invalid argument \"--N_iter\" for -evaluate");
					}
				}
				else if( (string(argv[carg_i]) == "--infer_only") or (string(argv[carg_i]) == "--not_infer")){
					if(string(argv[carg_i]) == "--infer_only"){
						infer_only = true;
					}
					else{
						no_infer = true;
					}

					//Now read the event nicknames and append to the list
					//Nicknames are added to the same list for infer_only or not_infer since they are exclusive
					//An exception will be raised later if both were given
					while( (carg_i+1<argc)
							and string(argv[carg_i+1]).size()>=1
							and string(argv[carg_i+1]).substr(0,1)!="-"){
						++carg_i;
						infer_restrict_nicknames.emplace(argv[carg_i]);
					}
				}
				else if(string(argv[carg_i]) == "--fix_err"){
					fix_err_rate = true;
				}


				else{
					return terminate_IGoR_with_error_message("Unknown argument \""+string(argv[carg_i])+"\" to specify inference/evaluate parameters");
				}
			}
		}



		else if(string(argv[carg_i]) == "-chain"){
			//Provide a boolean for the choice of a chain type (thus a set of genomic templates and model)
			chain_provided = true;
			++carg_i;
			if( (string(argv[carg_i]) == "alpha")
					or (string(argv[carg_i]) == "beta")
					or (string(argv[carg_i]) == "light")
					or (string(argv[carg_i]) == "heavy_naive")
					or (string(argv[carg_i]) == "heavy_memory")){
				chain_arg_str = string(argv[carg_i]);
				clog<<"Chain parameter set to: "<<chain_arg_str<<endl;
			}
			else{
				return terminate_IGoR_with_error_message("Unknown argument \""+string(argv[carg_i])+"\" to specify the chain choice!\n Supported arguments are: alpha, beta, heavy_naive , heavy_memory , light");
			}
		}

		else if(string(argv[carg_i]) == "-species"){
			//TODO add a check on the existence of the species
			species_provided = true;
			++carg_i;
			species_str = string(argv[carg_i]);
			clog<<"Species parameter set to: "<<species_str<<endl;
		}

		/*
		 * Set the working directory
		 * /!\ Needs to be set before all the counters !  /!\
		 */
		else if(string(argv[carg_i]) == "-set_wd"){
			wd = true;
			++carg_i;
			cl_path = string(string(argv[carg_i]));

			//Append a "/" if there is not one at the end of the directory path
			if (cl_path[cl_path.size()-1] != '/'){
				cl_path+="/";
			}

			if( not cl_counters_list.empty()){
				return terminate_IGoR_with_error_message("Working directory needs to be set before declaring the counters, please re-order the arguments");
			}
		}

		/*
		 * Output arguments parsing
		 */
		else if(string(argv[carg_i]) == "-output"){
			while( (carg_i+1<argc)
					and (string(argv[carg_i+1]).size()>2)
					and string(argv[carg_i+1]).substr(0,2)=="--"){

				++carg_i;
				/*
				 * TODO For now forget about outputing for every sequences / every iterations (more command line parameters to code)
				 */
				if(string(argv[carg_i]) == "--Pgen"){
					shared_ptr<Counter> pgen_counter_ptr(new Pgen_counter (cl_path + "output/"));
					cl_counters_list.emplace(cl_counters_list.size(),pgen_counter_ptr);
				}
				else if(string(argv[carg_i]) == "--scenarios"){
					int n_record_scenarios;
					++carg_i;
					try{
						n_record_scenarios = stoi(string(argv[carg_i]));
					}
					catch(exception& e){
						return terminate_IGoR_with_error_message("Expected the number of scenarios to be recorded by the best scenario counter, received: \"" + string(argv[carg_i]) + "\"");
					}

					if(n_record_scenarios<=0){
						return terminate_IGoR_with_error_message("Number of scenarios to be recorded must be greater than zero");
					}

					shared_ptr<Counter>best_sc_ptr(new Best_scenarios_counter(n_record_scenarios , cl_path + "output/" ,true));
					cl_counters_list.emplace(cl_counters_list.size(),best_sc_ptr);
				}
				else if(string(argv[carg_i]) == "--coverage"){
					Gene_class chosen_gc;
					++carg_i;
					try{
						chosen_gc = str2GeneClass(string(argv[carg_i]));
					}
					catch(exception& e){
						return terminate_IGoR_with_error_message("Unknown argument \""+string(argv[carg_i])+"\" to specify coverage target!\n Supported arguments are: V_gene, VD_genes, D_gene, DJ_gene, VJ_gene, J_gene, VDJ_genes");
					}
					shared_ptr<Counter> coverage_counter_ptr(new Coverage_err_counter(cl_path + "output/",chosen_gc,1,false,true));
					cl_counters_list.emplace(cl_counters_list.size(),coverage_counter_ptr);
				}
			}
		}

		/*
		 * Sequence generation arguments parsing
		 */
		else if(string(argv[carg_i]) == "-generate"){
			generate = true;
			++carg_i;

			//Number of sequences to generate must be given before optionnal arguments
			try{
				generate_n_seq = stoi(string(argv[carg_i]));
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Expected the number of sequences to generate, received: \"" + string(argv[carg_i]) + "\"");
			}

			while( (carg_i+1<argc)
					and (string(argv[carg_i+1]).size()>2)
					and string(argv[carg_i+1]).substr(0,2) == "--"){

				++carg_i;

				if(string(argv[carg_i]) == "--noerr"){
					generate_werr = false;
				}
				else if(string(argv[carg_i]) == "--CDR3"){
					gen_output_CDR3_data = true;
				}
				else if(string(argv[carg_i]) == "--name"){
					++carg_i;
					gen_filename_prefix = string(argv[carg_i])+"_";
				}
				else if(string(argv[carg_i]) == "--seed"){

					++carg_i;
					try{
						double tmp_seed = stod(string(argv[carg_i]));
						//Make sure the passed value is a positive integer
						if(tmp_seed<0
								or ((int) tmp_seed - tmp_seed)!=0.0){
							return terminate_IGoR_with_error_message("The random generator seed after \"--seed\" must be a positive integer");
						}
					}
					catch(exception& e){
						return terminate_IGoR_with_error_message("Expected a positive integer after \"--seed\" for sequence generation's seed, received: \"" + string(argv[carg_i]) + "\"");
					}
					gen_random_engine_seed = stoi(string(argv[carg_i]));
				}
				else{
					return terminate_IGoR_with_error_message("Unknown argument \""+string(argv[carg_i])+"\" to specify sequence generation parameters");
				}
			}

		}

		/*
		 * Input sequences argument parsing
		 */
		else if(string(argv[carg_i]) == "-read_seqs"){
			read_seqs = true;
			++carg_i;
			input_seqs_file = string(argv[carg_i]);
			//Get the extension index
			size_t extension_index = input_seqs_file.rfind(".");
			if(extension_index!=string::npos){
				string tmp_str = input_seqs_file.substr(extension_index , string::npos );
				transform(tmp_str.begin(),tmp_str.end(),tmp_str.begin(),::tolower);
				if( tmp_str == ".fasta" ){
					seqs_fileformat = FASTA_f;
					clog<<"FASTA extension detected for the input sequence file"<<endl;
				}
				else if(tmp_str == ".csv"){
					seqs_fileformat = CSV_f;
					clog<<"CSV extension detected for the input sequence file"<<endl;
				}
				else if(tmp_str == ".txt"){
					seqs_fileformat = TXT_f;
					clog<<"TXT extension detected for the input sequence file"<<endl;
				}
//				else if(tmp_str == ".fastq"){
//					seqs_fileformat = FASTQ_f;
//					clog<<"FASTQ extension detected for the input sequence file"<<endl;
//					clog<<"Parsing fastq to fasta "<<endl;
//					string basefilename = input_seqs_file.substr(0, extension_index);
//					try{
//					  system("fastaq_to_fasta -i "+input_seqs_file+" -o "+basefilename+".fasta");
//					}catch(exception& e){
//					  cout << e.what() <<endl;
//						return 1;
//					}
//					
//					clog<<"Rerun igor using the new fasta file createded :"<<basefilename <<".fasta"<<endl;
//				}
				else{
					return terminate_IGoR_with_error_message("Unknown file extension \"" + tmp_str + "\" for input sequences file! ");
				}
			}
			else{
				clog<<"No extension detected for the input sequence file assuming a text file without header"<<endl;
				seqs_fileformat = TXT_f;
			}
		}

		//If the argument does not correspond to any previous section throw an exception
		else{
			return terminate_IGoR_with_error_message("Unknown IGoR command line argument \""+string(argv[carg_i])+"\" ");
		}

		//Read the next command line argument
		++carg_i;
	}

	//Make sure the working directory is set somewhere before performing any action
	if(not wd){
		cl_path = "/tmp/";
	}
	clog<<"Working directory set to: \""+cl_path+"\""<<endl;

	//Check that both species and chain have been provided
	if(chain_provided xor species_provided){
		forward_list<string> error_messages;
		if(chain_provided){
			error_messages.emplace_front("Only chain argument was provided by the user");
		}
		else{
			error_messages.emplace_front("Only species argument was provided by the user");
		}
		error_messages.emplace_front("Both species and chain must be provided when using a predefined model!");
		return terminate_IGoR_with_error_message(error_messages);
	}


	if(chain_provided){
		if(chain_arg_str == "alpha"){
			has_D = false;
			chain_path_str = "tcr_alpha";
			try{
				v_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRA V genomic templates.",  e);
			}

			try{
				j_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRA J genomic templates.",  e);
			}
		}
		else if(chain_arg_str == "beta"){
			has_D = true;
			chain_path_str = "tcr_beta";
			try{
				v_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRB V genomic templates.",  e);
			}

			try{
				d_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicDs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRB D genomic templates.",  e);
			}

			try{
				j_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading TRB J genomic templates.",  e);
			}

		}
		else if(chain_arg_str == "light"){
			forward_list<string> error_messages;
			error_messages.emplace_front("If you wish to use IGoR on light chains please contact us so we can work on incorporating a light chain model to IGoR.");
			error_messages.emplace_front("Support for light chains in command line is not ready yet due to the lack of genomic templates and suitable model.");
			error_messages.emplace_front("Light chains support does not exist yet for command line!");
			return terminate_IGoR_with_error_message(error_messages);
		}
		else if( (chain_arg_str == "heavy_naive") or (chain_arg_str == "heavy_memory") ){
			has_D = true;
			chain_path_str = "bcr_heavy";
			try{
				v_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading IGH V genomic templates.",  e);
			}

			try{
				d_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicDs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading IGH D genomic templates.",  e);
			}

			try{
				j_genomic = read_genomic_fasta(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading IGH J genomic templates.",  e);
			}
			if( chain_arg_str == "heavy_naive" ){
				//Use a single error rate
			}
			else{
				//Memory

				//TODO infer only \mu for the hypermutation model
			}
		}
		//Read CDR3 anchors(cystein, tryptophan/phenylalanin indices)
		try{
			v_CDR3_anchors = read_gene_anchors_csv(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/V_gene_CDR3_anchors.csv");
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading V CDR3 anchors.",  e);
		}

		try{
			j_CDR3_anchors = read_gene_anchors_csv(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/J_gene_CDR3_anchors.csv");
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading J CDR3 anchors.",  e);
		}
	}

	//Read custom genomic templates if some custom ones were specified
	if(custom_v){
		try{
			v_genomic = read_genomic_fasta(custom_v_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom V genomic templates.",  e);
		}
	}
	if(custom_d){
		has_D = true;
		try{
			d_genomic = read_genomic_fasta(custom_d_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom D genomic templates.",  e);
		}
	}
	if(custom_j){
		try{
			j_genomic = read_genomic_fasta(custom_j_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom J genomic templates.",  e);
		}
	}

	//Read custom CDR3
	if(custom_v_anchors){
		try{
			v_CDR3_anchors = read_gene_anchors_csv(custom_v_anchors_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom V CDR3 anchors",  e);
		}
	}
	if(custom_j_anchors){
		try{
			j_CDR3_anchors = read_gene_anchors_csv(custom_j_anchors_path);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading user's custom J CDR3 anchors",  e);
		}
	}

	//Make sure passed arguments are unambiguous
	if(custom_cl_parms and load_last_inferred_parms){
		return terminate_IGoR_with_error_message("Setting a custom model and loading the last inferred model in the same command is ambiguous!");
	}

	//Load last inferred model
	if(load_last_inferred_parms){
		clog<<"Loading last inferred model..."<<endl;
		try{
			cl_model_parms.read_model_parms(cl_path +  batchname + "inference/final_parms.txt");
			cl_model_marginals = Model_marginals(cl_model_parms);
			cl_model_marginals.txt2marginals(cl_path +  batchname + "inference/final_marginals.txt",cl_model_parms);

			//Check if the model contains a D gene event in order to load the alignments
			auto events_map = cl_model_parms.get_events_map();
			if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side))>0){
				has_D = true;
			}
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while loading last inferred model, please check that the model exists",e);
		}
	}

	/*
	 * Read supplied model parms and marginals
	 */
	if( ((not custom_cl_parms) and (not load_last_inferred_parms))
			and (infer or evaluate or generate)){
		clog<<"Read some model parms"<<endl;
		try{
			cl_model_parms.read_model_parms(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/models/model_parms.txt");
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading IGoR's model parameters.",e);
		}
		cl_model_marginals = Model_marginals(cl_model_parms);
		try{
			cl_model_marginals.txt2marginals(string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/models/model_marginals.txt",cl_model_parms);
		}
		catch(exception& e){
			return terminate_IGoR_with_error_message("Exception caught while reading IGoR's model marginals.",e);
		}
	}


	/*
	 * If some custom genomic templates were supplied, two possible cases here:
	 * - if any supplied genomic template is absent from the model, or if its actual sequence is different the marginals will be re-initialized
	 * - if all the supplied templates were already contained in the model the missing one will be set to 0 probability,
	 * 	 the others will keep their probability ratio.
	 *
	 * This will be executed whether using a supplied model, a custom model or the last inferred one.
	 */
	if((infer or evaluate or generate)){
		bool any_custom_gene = false;
		unordered_map<tuple<Event_type,Gene_class,Seq_side>,shared_ptr<Rec_Event>> tmp_events_map = cl_model_parms.get_events_map();
		if(custom_v){
			shared_ptr<Rec_Event> v_choice = tmp_events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side));
			shared_ptr<Gene_choice> v_choice_gc = dynamic_pointer_cast<Gene_choice>(v_choice);
			bool any_genomic_difference = false;
			unordered_map<string , Event_realization> realization_map_copy = v_choice_gc->get_realizations_map();
			/*
			 * We loop over provided genomic templates and check if they are contained in the current model
			 */
			for(pair<string,string> genomic_template : v_genomic){
				if(realization_map_copy.count(genomic_template.first)>0){
					Event_realization& ev_real = realization_map_copy.at(genomic_template.first);
					if(ev_real.value_str == genomic_template.second){
						//Remove the genomic template from the unseen list
						//FIXME this will fail if the same genomic template has been read several times
						realization_map_copy.erase(genomic_template.first);
					}
					else{
						any_genomic_difference = true;
					}
				}
				else{
					any_genomic_difference = true;
				}
				if(any_genomic_difference){
					break;
				}
			}
			if(any_genomic_difference){
				//If any difference is detected reset the marginals
				Rec_Event_name former_name = v_choice_gc->get_name();
				v_choice_gc->set_genomic_templates(v_genomic);
				cl_model_parms.update_edge_event_name(former_name,v_choice_gc->get_name());
				any_custom_gene = true;
			}
			else{
				//Else we set to 0 probability all the ones that were not found in the list of genomic templates
				for(pair<string,Event_realization> ev_real : realization_map_copy){
					cl_model_marginals.set_realization_proba(ev_real.first,v_choice,0.0,cl_model_parms);
				}
			}
		}
		if(has_D and custom_d){
			shared_ptr<Rec_Event> d_choice = tmp_events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side));
			shared_ptr<Gene_choice> d_choice_gc = dynamic_pointer_cast<Gene_choice>(d_choice);
			bool any_genomic_difference = false;
			unordered_map<string , Event_realization> realization_map_copy = d_choice_gc->get_realizations_map();
			/*
			 * We loop over provided genomic templates and check if they are contained in the current model
			 */
			for(pair<string,string> genomic_template : d_genomic){
				if(realization_map_copy.count(genomic_template.first)>0){
					Event_realization& ev_real = realization_map_copy.at(genomic_template.first);
					if(ev_real.value_str == genomic_template.second){
						//Remove the genomic template from the unseen list
						//FIXME this will fail if the same genomic template has been read several times
						realization_map_copy.erase(genomic_template.first);
					}
					else{
						any_genomic_difference = true;
					}
				}
				else{
					any_genomic_difference = true;
				}
				if(any_genomic_difference){
					break;
				}
			}
			if(any_genomic_difference){
				//If any difference is detected reset the marginals
				Rec_Event_name former_name = d_choice_gc->get_name();
				d_choice_gc->set_genomic_templates(d_genomic);
				cl_model_parms.update_edge_event_name(former_name,d_choice_gc->get_name());
				any_custom_gene = true;
			}
			else{
				//Else we set to 0 probability all the ones that were not found in the list of genomic templates
				for(pair<string,Event_realization> ev_real : realization_map_copy){
					cl_model_marginals.set_realization_proba(ev_real.first,d_choice,0.0,cl_model_parms);
				}
			}
		}
		if(custom_j){
			shared_ptr<Rec_Event> j_choice = tmp_events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side));
			shared_ptr<Gene_choice> j_choice_gc = dynamic_pointer_cast<Gene_choice>(j_choice);
			bool any_genomic_difference = false;
			unordered_map<string , Event_realization> realization_map_copy = j_choice_gc->get_realizations_map();
			/*
			 * We loop over provided genomic templates and check if they are contained in the current model
			 */
			for(pair<string,string> genomic_template : j_genomic){
				if(realization_map_copy.count(genomic_template.first)>0){
					Event_realization& ev_real = realization_map_copy.at(genomic_template.first);
					if(ev_real.value_str == genomic_template.second){
						//Remove the genomic template from the unseen list
						//FIXME this will fail if the same genomic template has been read several times
						realization_map_copy.erase(genomic_template.first);
					}
					else{
						any_genomic_difference = true;
					}
				}
				else{
					any_genomic_difference = true;
				}
				if(any_genomic_difference){
					break;
				}
			}
			if(any_genomic_difference){
				//If any difference is detected reset the marginals
				Rec_Event_name former_name = j_choice_gc->get_name();
				j_choice_gc->set_genomic_templates(j_genomic);
				cl_model_parms.update_edge_event_name(former_name,j_choice_gc->get_name());
				any_custom_gene = true;
			}
			else{
				//Else we set to 0 probability all the ones that were not found in the list of genomic templates
				for(pair<string,Event_realization> ev_real : realization_map_copy){
					cl_model_marginals.set_realization_proba(ev_real.first,j_choice,0.0,cl_model_parms);
				}
			}
		}

		if(any_custom_gene){
			/*
			 * If some custom genomic templates were provided we have replaced the genomic templates contained in the model
			 * Thus other components (e.g deletion profiles) might not match anymore and re inferring a model is necessary
			 * Marginals are initialized with a uniform distribution
			 */
			cl_model_marginals = Model_marginals(cl_model_parms);
			cl_model_marginals.uniform_initialize(cl_model_parms);
			clog<<"Not all custom genomic templates were found in the loaded Model parameters, Model marginals are thus reinitialized to a uniform distribution!"<<endl;
		}

	}

	/*
	 * Once model parms have been read fix the events requested
	 */
	if(infer_only and no_infer){
		return terminate_IGoR_with_error_message("Cannot use both \"--infer_only\" and \"--not_infer\" since they are somewhat redundant");
	}
	if((infer_only or no_infer)
		and (infer or evaluate or generate)){
		if(infer_only){
			//Loop over events and fix all but the ones given in the list
			list<shared_ptr<Rec_Event>> events_list = cl_model_parms.get_event_list();
			for(list<shared_ptr<Rec_Event>>::iterator iter = events_list.begin() ; iter!=events_list.end() ; ++iter){
				if(infer_restrict_nicknames.count((*iter)->get_nickname()) >0){
					(*iter)->fix(false); //Technically not useful since set to false by default, just a safety
				}
				else{
					(*iter)->fix(true);
				}
			}
			//TODO add a check to see whether all nicknames exist
		}
		else{
			//Fix all events provided
			for(set<string>::const_iterator iter = infer_restrict_nicknames.begin() ; iter!=infer_restrict_nicknames.end() ; ++iter){
				shared_ptr<Rec_Event>event_ptr = cl_model_parms.get_event_pointer((*iter),true);
				event_ptr->fix(true);
			}
		}
	}

	/*
	 * Fix the error rate if requested
	 */
	if(fix_err_rate
		and (infer or evaluate or generate)){
		cl_model_parms.get_err_rate_p()->update_value(false);
	}

	//If a batchname is defined edit the output dir name
	//This might be overcomplicated for no reason but should be easy to maintain
	for(pair<size_t,shared_ptr<Counter>> counter_int_pair_ptr : cl_counters_list){
		string former_path = counter_int_pair_ptr.second->get_path_to_files();
		size_t output_str_index = former_path.rfind("/output/");
		string new_path = former_path.substr(0,output_str_index+1) + batchname + "output" + former_path.substr(output_str_index+7,string::npos); //7 is the length of /output
		counter_int_pair_ptr.second->set_path_to_files(new_path);
	}

	//Output warnings on the use of the subsample command as it may have different effects depending on the command used
	if(subsample_seqs){
		if(read_seqs){
			clog<<"Subsampling "<<n_subsample_seqs<<" sequences from the input sequence file, the resulting indexed sequence will be a subsample."<<endl;
		}
		else if(align){
			if((infer or evaluate)){
				clog<<"Subsampling "<<n_subsample_seqs<<" sequences for alignments and evaluation/inference"<<endl;
			}
			else{
				clog<<"Subsampling "<<n_subsample_seqs<<" sequences for alignments without an -evaluate or -infer command supplied."<<endl;
				clog<<"/!\\ Because -subsample N makes a random sample do not further re-use the command with -evaluate or -infer /!\\ "<<endl;
			}
		}
		else if((infer or evaluate)){
			clog<<"Subsampling "<<n_subsample_seqs<<" sequences for evaluation/inference without an -align command supplied."<<endl;
			clog<<"/!\\ Because -subsample N makes a random sample make sure you have not used this command during alignment, otherwise the resulting sample will be the intersection of the two random samples! The resulting sample size cannot be guaranteed /!\\ "<<endl;
		}
	}




	if(run_demo){

		/*Run this sample demo code
		 *
		 * Outline:
		 *
		 * Read TCRb genomic templates
		 *
		 * Align the sequences contained in the /demo/murugan_naive1_noncoding_demo_seqs.txt file to those templates
		 *
		 * Create a TCRb model, a simple error rate and the corresponding marginals
		 *
		 * Show reads and write functions for the model and marginals
		 *
		 * Infer a model from the sequences (perform 10 iterations of EM)
		 *
		 * Generate sequences from the obtained model
		 *
		 */

		//Path to the working folder
/*		string path (argv[1]);
		if (path[path.size()-1] != '/'){
			path+="/";
		}*/

		clog<<"Reading genomic templates"<<endl;

		vector<pair<string,string>> v_genomic = read_genomic_fasta( string(IGOR_DATA_DIR) + "/demo/genomicVs_with_primers.fasta");

		vector<pair<string,string>> d_genomic = read_genomic_fasta( string(IGOR_DATA_DIR) + "/demo/genomicDs.fasta");

		vector<pair<string,string>> j_genomic = read_genomic_fasta( string(IGOR_DATA_DIR) + "/demo/genomicJs_all_curated.fasta");

		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = { // A,C,G,T,R,Y,K,M,S,W,B,D,H,V,N
		        5,-14,-14,-14,-14,2,-14,2,2,-14,-14,1,1,1,0,
		        -14,5,-14,-14,-14,2,2,-14,-14,2,1,-14,1,1,0,
		        -14,-14,5,-14,2,-14,2,-14,2,-14,1,1,-14,1,0,
		        -14,-14,-14,5,2,-14,-14,2,-14,2,1,1,1,-14,0,
		        -14,-14,2,2,1.5,-14,-12,-12,-12,-12,1,1,-13,-13,0,
		        2,2,-14,-14,-14,1.5,-12,-12,-12,-12,-13,-13,1,1,0,
		        -14,2,2,-14,-12,-12,1.5,-14,-12,-12,1,-13,-13,1,0,
		        2,-14,-14,2,-12,-12,-14,1.5,-12,-12,-13,1,1,-13,0,
		        2,-14,2,-14,-12,-12,-12,-12,1.5,-14,-13,1,-13,1,0,
		        -14,2,-14,2,-12,-12,-12,-12,-14,1.5,1,-13,1,-13,0,
		        -14,1,1,1,1,-13,1,-13,-13,1,0.5,-12,-12,-12,0,
		        1,-14,1,1,1,-13,-13,1,1,-13,-12,0.5,-12,-12,0,
		        1,1,-14,1,-13,1,-13,1,-13,1,-12,-12,0.5,-12,0,
		        1,1,1,-14,-13,1,1,-13,1,-13,-12,-12,-12,0.5,0,
		        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		Matrix<double> nuc44_sub_matrix(15,15,nuc44_vect);

		//Instantiate aligner with substitution matrix declared above and gap penalty of 50
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
		d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		clog<<"Reading sequences and aligning"<<endl;
		typedef std::chrono::system_clock myclock;
		myclock::time_point begin_time, end_time;

		begin_time = myclock::now();


		vector<pair<const int, const string>> indexed_seqlist = read_txt( string(IGOR_DATA_DIR) + "/demo/murugan_naive1_noncoding_demo_seqs.txt" ); //Could also read a FASTA file <code>read_fasta()<\code> or indexed sequences <code>read_indexed_seq_csv()<\code>

		cl_path+="igor_demo/";
		system(&("mkdir " + cl_path )[0]);

		v_aligner.align_seqs( string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_V.csv"),indexed_seqlist,50,true,INT16_MIN,-155);
		//v_aligner.write_alignments_seq_csv(path + string("alignments_V.csv") , v_alignments);

		d_aligner.align_seqs(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_D.csv"),indexed_seqlist,0,false);
		//d_aligner.write_alignments_seq_csv(path + string("alignments_D.csv") , d_alignments);

		j_aligner.align_seqs(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_J.csv"),indexed_seqlist,10,true,42,48);

		end_time= myclock::now();
		chrono::duration<double> elapsed = end_time - begin_time;
		clog<<"Alignments procedure lasted: "<<elapsed.count()<<" seconds"<<endl;
		clog<<"for "<<indexed_seqlist.size()<<" TCRb sequences of 60bp(from murugan and al), against ";
		clog<<v_genomic.size()<<" Vs,"<<d_genomic.size()<<" Ds, and "<<j_genomic.size()<<" Js full sequences"<<endl;

		//unordered_map<int,forward_list<Alignment_data>> j_alignments = j_aligner.align_seqs(indexed_seqlist,10,true,42,48);
		//j_aligner.write_alignments_seq_csv(path + string("alignments_J.csv") , j_alignments);


		write_indexed_seq_csv(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_indexed_seq.csv") , indexed_seqlist);



		clog<<"Construct the model"<<endl;
		//Construct a TCRb model
		Gene_choice v_choice(V_gene,v_genomic);
		v_choice.set_nickname("v_choice");
		v_choice.set_priority(7);
		Gene_choice d_choice(D_gene,d_genomic);
		d_choice.set_nickname("d_gene");
		d_choice.set_priority(6);
		Gene_choice j_choice(J_gene,j_genomic);
		j_choice.set_nickname("j_choice");
		j_choice.set_priority(7);

		Deletion v_3_del(V_gene,Three_prime,make_pair(-4,16));//16
		v_3_del.set_nickname("v_3_del");
		v_3_del.set_priority(5);
		Deletion d_5_del(D_gene,Five_prime,make_pair(-4,16));
		d_5_del.set_nickname("d_5_del");
		d_5_del.set_priority(5);
		Deletion d_3_del(D_gene,Three_prime,make_pair(-4,16));
		d_3_del.set_nickname("d_3_del");
		d_3_del.set_priority(5);
		Deletion j_5_del(J_gene,Five_prime,make_pair(-4,18));
		j_5_del.set_nickname("j_5_del");
		j_5_del.set_priority(5);

		Insertion vd_ins(VD_genes,make_pair(0,30));
		vd_ins.set_nickname("vd_ins");
		vd_ins.set_priority(4);
		Insertion dj_ins(DJ_genes,make_pair(0,30));
		dj_ins.set_nickname("dj_ins");
		dj_ins.set_priority(2);

		Dinucl_markov markov_model_vd(VD_genes);
		markov_model_vd.set_nickname("vd_dinucl");
		markov_model_vd.set_priority(3);

		Dinucl_markov markov_model_dj(DJ_genes);
		markov_model_dj.set_nickname("dj_dinucl");
		markov_model_dj.set_priority(1);


		Model_Parms parms;

		//Add nodes to the graph
		parms.add_event(&v_choice);
		parms.add_event(&d_choice);
		parms.add_event(&j_choice);

		parms.add_event(&v_3_del);
		parms.add_event(&d_3_del);
		parms.add_event(&d_5_del);
		parms.add_event(&j_5_del);

		parms.add_event(&vd_ins);
		parms.add_event(&dj_ins);

		parms.add_event(&markov_model_vd);
		parms.add_event(&markov_model_dj);


		//Add correlations
		parms.add_edge(&v_choice,&v_3_del);
		parms.add_edge(&j_choice,&j_5_del);
		parms.add_edge(&d_choice,&d_3_del);
		parms.add_edge(&d_choice,&d_5_del);
		parms.add_edge(&d_5_del,&d_3_del);
		parms.add_edge(&j_choice,&d_choice);


		//Create the corresponding marginals
		Model_marginals model_marginals(parms);
		model_marginals.uniform_initialize(parms); //Can also start with a random prior using random_initialize()

		//Instantiate an error rate
		Single_error_rate error_rate(0.001);

		parms.set_error_ratep(&error_rate);

		clog<<"Write and read back the model"<<endl;
		//Write the model_parms into a file
		parms.write_model_parms(string(cl_path + "/demo_write_model_parms.txt"));

		//Write the marginals into a file
		model_marginals.write2txt(string(cl_path + "/demo_write_model_marginals.txt"),parms);

		//Read a model and marginal pair
		Model_Parms read_model_parms;
		read_model_parms.read_model_parms(string(cl_path + "/demo_write_model_parms.txt"));
		Model_marginals read_model_marginals(read_model_parms);
		read_model_marginals.txt2marginals(string(cl_path + "/demo_write_model_marginals.txt"),read_model_parms);

		//Instantiate a Counter
		map<size_t,shared_ptr<Counter>> counters_list;
		 //Collect gene coverage and errors
		shared_ptr<Counter> coverage_counter_ptr(new Coverage_err_counter(cl_path + "/run_demo/",VJ_genes,1,false,false));
		counters_list.emplace(0,coverage_counter_ptr);

		 //Collect 10 best scenarios per sequence during the last iteration
		shared_ptr<Counter>best_sc_ptr(new Best_scenarios_counter(10 , cl_path + "/run_demo/" ,true ));
		counters_list.emplace(1,best_sc_ptr);

		 //Collect sequence generation probability during last iteration
		shared_ptr<Counter> pgen_counter_ptr(new Pgen_counter (cl_path + "/run_demo/"));
		counters_list.emplace(2,pgen_counter_ptr);

		shared_ptr<Counter> errors_counter(new Errors_counter (10,string(cl_path + "/run_demo/")));
		counters_list.emplace(3,errors_counter);

		//Instantiate the high level GenModel class
		//This class allows to make most useful high level operations(model inference/Pgen computation , sequence generation)
		GenModel gen_model(read_model_parms,read_model_marginals,counters_list);


		//Inferring a model

		//Read alignments
		//vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path+ string(argv[2]) + string("indexed_seq.csv"));
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
		sorted_alignments = read_alignments_seq_csv_score_range(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15
		sorted_alignments = read_alignments_seq_csv_score_range(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<tuple<int,string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

		//Infer the model
		clog<<"Infer model"<<endl;

		begin_time = myclock::now();
		system(&("mkdir " + cl_path + "run_demo")[0]);
		gen_model.infer_model(sorted_alignments_vec , 4 , string(cl_path + "/run_demo/") , true ,1e-35,0.0001);

		end_time= myclock::now();
		elapsed = end_time - begin_time;
		clog<<"Model inference procedure lasted: "<<elapsed.count()<<" seconds"<<endl;


		//Generate sequences
		clog<<"Generate sequences"<<endl;
		auto generated_seq =  gen_model.generate_sequences(100,false);//Without errors
		gen_model.write_seq_real2txt(string(cl_path + "/generated_seqs_indexed_demo.csv"), string(cl_path + "/generated_seqs_realizations_demo.csv") , generated_seq);//Member function will be changed

	}

	else if (not custom){
		//Execute code dictated by command line arguments
		if(read_seqs){
			vector<pair<const int, const string>> indexed_seqlist;
			switch(seqs_fileformat){
			case FASTA_f:
				indexed_seqlist = read_fasta(input_seqs_file);
				break;
			case CSV_f:
				indexed_seqlist = read_indexed_csv(input_seqs_file);
				break;
			case TXT_f:
				indexed_seqlist = read_txt(input_seqs_file);
				break;
			default:
				return terminate_IGoR_with_error_message("Unknown file format for input seqs file!");
			}

			//create the directory
			system(&("mkdir " + cl_path + "aligns")[0]);

			if(subsample_seqs){
				try{
					indexed_seqlist = sample_indexed_seq(indexed_seqlist,n_subsample_seqs);
				}
				catch(exception& e){
					forward_list<string> error_messages;
					return terminate_IGoR_with_error_message("Exception caught trying to subsample input files sequences to indexed sequences:",e);
				}
			}

			write_indexed_seq_csv(cl_path + "aligns/" + batchname + "indexed_sequences.csv",indexed_seqlist);
		}

		if(align){
			vector<pair<const int, const string>> indexed_seqlist;
			try{
				indexed_seqlist = read_indexed_csv(cl_path + "aligns/" + batchname + "indexed_sequences.csv");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading indexed sequences file sequences before alignment. Make sure indexed sequence file has previously been created using \"-read_seqs\" with similar path parameters (working directory, batchname, ...)",e);
			}

			if(subsample_seqs and not read_seqs){
				try{
					indexed_seqlist = sample_indexed_seq(indexed_seqlist,n_subsample_seqs);
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught trying to subsample input files sequences before alignment:",e);
				}
			}

			if(align_v){
				//Performs V alignments
				Aligner v_aligner = Aligner(v_subst_matrix , v_gap_penalty , V_gene);
				v_aligner.set_genomic_sequences(v_genomic);
				try{
					if (not align_data_is_CDR3){
						clog<<"Performing V alignments...."<<endl;
						if(v_template_bounds_map.empty()){
							v_aligner.align_seqs(cl_path + "aligns/" +  batchname + v_align_filename , indexed_seqlist , v_align_thresh_value , v_best_align_only , v_best_gene_only , v_left_offset_bound , v_right_offset_bound,v_reversed_offsets);
						}
						else{
							// Check if all genes have a specific template entry
							// If not add the min max offsets bounds as entries
							list<string> unknown_template_specific_offsets;
							for(pair<string,string> v_template: v_genomic){
								if(v_template_bounds_map.count(v_template.first)==0){
									v_template_bounds_map.emplace(v_template.first,make_pair(v_left_offset_bound,v_right_offset_bound));
									unknown_template_specific_offsets.push_front(v_template.first);
								}
							}
							if(not unknown_template_specific_offsets.empty()){
								clog<<"Template specific offsets could not be found for the following "<<unknown_template_specific_offsets.size()<<" genes: ";
								for(string unknown_gene: unknown_template_specific_offsets) clog<<"\""<<unknown_gene<<"\" ,";
								clog<<endl<<"For these genes provided/default values for V gene min and max offset bounds ("<<v_left_offset_bound<<"/"<<v_right_offset_bound<<") been set as genomic offset bounds."<<endl;
							}
							v_aligner.align_seqs(cl_path + "aligns/" +  batchname + v_align_filename , indexed_seqlist , v_align_thresh_value , v_best_align_only , v_best_gene_only , v_template_bounds_map,v_reversed_offsets);
						}
					}
					else{ //Provided sequences are ntCDR3s, alignment offsets are based on provided CDR3 gene anchors.
						clog<<"Performing CDR3s V alignments...."<<endl;
						unordered_map<string,pair<int,int>> v_genomic_CDR3_offset_bounds;
						list<string> unknown_gene_anchors;
						int min_offset = INT32_MAX;
						int max_offset = INT32_MIN;
						for(pair<string,string> v_template: v_genomic){
							if(v_CDR3_anchors.count(v_template.first)>0){
								int gene_offset = - v_CDR3_anchors.at(v_template.first);
								//Use a reversed offset and substract 2 in order to take into account the anchor's codon
								v_genomic_CDR3_offset_bounds.emplace(v_template.first, make_pair(gene_offset,gene_offset));
								//Update min/max offsets values
								if(min_offset>gene_offset) min_offset = gene_offset;
								else if(max_offset<gene_offset) max_offset = gene_offset;
							}
							else{
								if(v_template_bounds_map.count(v_template.first)){
									v_genomic_CDR3_offset_bounds.emplace(v_template.first, v_template_bounds_map.at(v_template.first));
								}
								else{
									v_genomic_CDR3_offset_bounds.emplace(v_template.first, make_pair(v_left_offset_bound,v_right_offset_bound));
								}
								unknown_gene_anchors.push_front(v_template.first);
								/*Note: there was a tradeoff between
								 * - putting everything between the same min and max (extracted from template specific offsets) => most efficient, likely to be correct but dangerous
								 * - put no constraint => safest in theory, in practice risky if only best alignments are accepted.
								 * In the end this solution is more elegant and encapsulate both.
								 */
							}
						}
						if(not unknown_gene_anchors.empty()){
							clog<<"Anchors indices could not be found for the following "<<unknown_gene_anchors.size()<<" genes: ";
							for(string unknown_gene: unknown_gene_anchors) clog<<"\""<<unknown_gene<<"\" ,";
							clog<<endl<<"For these genes provided/default values for V gene min and max offset bounds ("<<v_left_offset_bound<<"/"<<v_right_offset_bound<<") or, if provided, template specific ones have been set as genomic offset bounds."<<endl;
							clog<<"Hint: provided CDR3 anchors correspond to min/max offsets in ["<<min_offset<<":"<<max_offset<<"]."<<endl;
							clog<<"If you have provided V min/max offset values make sure they were NOT defined as reversed offsets."<<endl;
						}
						//Call the aligner module
						v_aligner.align_seqs(cl_path + "aligns/" +  batchname + v_align_filename , indexed_seqlist , v_align_thresh_value , v_best_align_only , v_best_gene_only , v_genomic_CDR3_offset_bounds,false);
					}

					// Call the aligner
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught upon aligning V genomic templates.",e);
				}
			}


			if(has_D and align_d){
				//Performs D alignments if the chain contains a D
				clog<<"Performing D alignments...."<<endl;
				Aligner d_aligner = Aligner(d_subst_matrix , d_gap_penalty , D_gene);
				d_aligner.set_genomic_sequences(d_genomic);
				try{
					if(d_template_bounds_map.empty()){
						d_aligner.align_seqs(cl_path + "aligns/" +  batchname + d_align_filename ,indexed_seqlist, d_align_thresh_value , d_best_align_only , d_best_gene_only , d_left_offset_bound , d_right_offset_bound, d_reversed_offsets);
					}
					else{
						// Check if all genes have a specific template entry
						// If not add the min max offsets bounds as entries
						list<string> unknown_template_specific_offsets;
						for(pair<string,string> d_template: d_genomic){
							if(d_template_bounds_map.count(d_template.first)==0){
								d_template_bounds_map.emplace(d_template.first,make_pair(d_left_offset_bound,d_right_offset_bound));
								unknown_template_specific_offsets.push_front(d_template.first);
							}
						}
						if(not unknown_template_specific_offsets.empty()){
							clog<<"Template specific offsets could not be found for the following "<<unknown_template_specific_offsets.size()<<" genes: ";
							for(string unknown_gene: unknown_template_specific_offsets) clog<<"\""<<unknown_gene<<"\" ,";
							clog<<endl<<"For these genes provided/default values for D gene min and max offset bounds ("<<v_left_offset_bound<<"/"<<v_right_offset_bound<<") been set as genomic offset bounds."<<endl;
						}
						d_aligner.align_seqs(cl_path + "aligns/" +  batchname + d_align_filename , indexed_seqlist , d_align_thresh_value , d_best_align_only , d_best_gene_only , d_template_bounds_map, d_reversed_offsets);
					}
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught upon aligning D genomic templates.",e);
				}
			}

			if(align_j){
				Aligner j_aligner (j_subst_matrix , j_gap_penalty , J_gene);
				j_aligner.set_genomic_sequences(j_genomic);
				try{
					if (not align_data_is_CDR3){
						clog<<"Performing J alignments...."<<endl;
						if(j_template_bounds_map.empty()){
							j_aligner.align_seqs(cl_path + "aligns/" +  batchname + j_align_filename , indexed_seqlist, j_align_thresh_value , j_best_align_only , j_best_gene_only , j_left_offset_bound , j_right_offset_bound, j_reversed_offsets);
						}
						else{
							// Check if all genes have a specific template entry
							// If not add the min max offsets bounds as entries
							list<string> unknown_template_specific_offsets;
							for(pair<string,string> j_template: j_genomic){
								if(j_template_bounds_map.count(j_template.first)==0){
									j_template_bounds_map.emplace(j_template.first,make_pair(j_left_offset_bound,j_right_offset_bound));
									unknown_template_specific_offsets.push_front(j_template.first);
								}
							}
							if(not unknown_template_specific_offsets.empty()){
								clog<<"Template specific offsets could not be found for the following "<<unknown_template_specific_offsets.size()<<" genes: ";
								for(string unknown_gene: unknown_template_specific_offsets) clog<<"\""<<unknown_gene<<"\" ,";
								clog<<endl<<"For these genes provided/default values for J gene min and max offset bounds ("<<v_left_offset_bound<<"/"<<v_right_offset_bound<<") been set as genomic offset bounds."<<endl;
							}
							j_aligner.align_seqs(cl_path + "aligns/" +  batchname + j_align_filename , indexed_seqlist , j_align_thresh_value , j_best_align_only , j_best_gene_only , j_template_bounds_map, j_reversed_offsets);
						}
					}
					else{ //Provided sequences are ntCDR3s, alignment offsets are based on provided CDR3 gene anchors.
						clog<<"Performing CDR3s J alignments...."<<endl;
						unordered_map<string,pair<int,int>> j_genomic_CDR3_offset_bounds;
						list<string> unknown_gene_anchors;
						int min_offset = INT32_MAX;
						int max_offset = INT32_MIN;
						for(pair<string,string> j_template: j_genomic){
							if(j_CDR3_anchors.count(j_template.first)>0){
								int gene_offset = - j_CDR3_anchors.at(j_template.first) -2;
								//Use a reversed offset and substract 2 in order to take into account the anchor's codon
								j_genomic_CDR3_offset_bounds.emplace(j_template.first, make_pair(gene_offset,gene_offset));
								//Update min/max offsets values
								if(min_offset>gene_offset) min_offset = gene_offset;
								else if(max_offset<gene_offset) max_offset = gene_offset;
							}
							else{
								if(j_template_bounds_map.count(j_template.first)){
									j_genomic_CDR3_offset_bounds.emplace(j_template.first, j_template_bounds_map.at(j_template.first));
								}
								else{
									j_genomic_CDR3_offset_bounds.emplace(j_template.first, make_pair(j_left_offset_bound,j_right_offset_bound));
								}
								unknown_gene_anchors.push_front(j_template.first);
								/*Note: there was a tradeoff between
								 * - putting everything between the same min and max (extracted from template specific offsets) => most efficient, likely to be correct but dangerous
								 * - put no constraint => safest in theory, in practice risky if only best alignments are accepted.
								 * In the end this solution is more elegant and encapsulate both.
								 */
							}
						}
						if(not unknown_gene_anchors.empty()){
							clog<<"Anchors indices could not be found for the following "<<unknown_gene_anchors.size()<<" genes: ";
							for(string unknown_gene: unknown_gene_anchors) clog<<"\""<<unknown_gene<<"\" ,";
							clog<<endl<<"For these genes provided/default values for J gene min and max offset bounds ("<<j_left_offset_bound<<"/"<<j_right_offset_bound<<") or, if provided, template specific ones have been set as genomic offset bounds."<<endl;;
							clog<<"Hint: provided CDR3 anchors correspond to min/max offsets in ["<<min_offset<<":"<<max_offset<<"]."<<endl;
							clog<<"If you have provided J min/max offset or gene specific offsets values make sure they were defined as reversed offsets."<<endl;
						}
						//Call the aligner module
						j_aligner.align_seqs(cl_path + "aligns/" +  batchname + j_align_filename , indexed_seqlist, j_align_thresh_value , j_best_align_only , j_best_gene_only , j_genomic_CDR3_offset_bounds,true);
					}
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught upon aligning J genomic templates.",e);
				}
			}
			
			//Get CDR3 from alignments.
			if (bFeature_CDR3){
				unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments;
				try{
					sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + v_align_filename, V_gene , 55 , false , indexed_seqlist  );
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught while reading V alignments before inference/evaluation. Make sure alignments were carried previously using \"-align --V\" or \"-align --all\" with similar path parameters (working directory, batchname, ...)",e);
				}

				try{
					sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + j_align_filename, J_gene , 10 , false , indexed_seqlist , sorted_alignments);
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught while reading J alignments before inference/evaluation. Make sure alignments were carried previously using \"-align --J\" or \"-align --all\" with similar path parameters (working directory, batchname, ...)",e);
				}
			
				string cl_path_ref_genome = string(IGOR_DATA_DIR) + "/models/"+species_str+"/"+chain_path_str+"/ref_genome/";
				string cl_path_aligns     = cl_path + "aligns/" +batchname;

				string flnV_CDR3_anchors = cl_path_ref_genome + "V_gene_CDR3_anchors.csv";
				string flnJ_CDR3_anchors = cl_path_ref_genome + "J_gene_CDR3_anchors.csv";


				string flnIndexedCDR3      = cl_path_aligns     + "indexed_CDR3.csv";
				ofstream ofileIndexedCDR3(flnIndexedCDR3);
				ofileIndexedCDR3 << "seq_index;v_anchor;j_anchor;CDR3nt;CDR3aa"<<endl;
			
				ExtractFeatures featureCDR3;
				featureCDR3.load_VJgenomicTemplates(v_genomic, j_genomic);
				featureCDR3.load_VJanchors(flnV_CDR3_anchors, flnJ_CDR3_anchors);
				featureCDR3.set_sorted_alignments(&sorted_alignments);
				
				// For each sequence get the CDR3
				for (auto seq_it = indexed_seqlist.begin(); seq_it != indexed_seqlist.end(); ++seq_it){
					CDR3SeqData cdr3InputSeq;
					int seq_index = (*seq_it).first;
					cdr3InputSeq = featureCDR3.generateCDR3( seq_index );
					ofileIndexedCDR3 << featureCDR3.generateCDR3_csv_line(cdr3InputSeq) << endl;
				}
				ofileIndexedCDR3.close();
			} // end extractCDR3
		}//end align

		if(infer xor evaluate){

			GenModel genmodel(cl_model_parms,cl_model_marginals,cl_counters_list);

			//Reading alignments
			vector<pair<const int, const string>> indexed_seqlist;
			try{
				indexed_seqlist = read_indexed_csv(cl_path + "aligns/" + batchname + "indexed_sequences.csv");
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading indexed sequences file sequences before inference/evaluation. Make sure indexed sequence file has previously been created using \"-read_seqs\" with similar path parameters (working directory, batchname, ...)",e);
			}

			if(subsample_seqs and not align){
				try{
					indexed_seqlist = sample_indexed_seq(indexed_seqlist,n_subsample_seqs);
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught trying to subsample indexed sequences before inference/evaluation:",e);
				}
			}
			unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments;
			try{
				sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + v_align_filename, V_gene , 55 , false , indexed_seqlist  );
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading V alignments before inference/evaluation. Make sure alignments were carried previously using \"-align --V\" or \"-align --all\" with similar path parameters (working directory, batchname, ...)",e);
			}

			if(has_D){
				try{
					sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + d_align_filename, D_gene , 35 , false , indexed_seqlist , sorted_alignments);
				}
				catch(exception& e){
					return terminate_IGoR_with_error_message("Exception caught while reading D alignments before inference/evaluation. Make sure alignments were carried previously using \"-align --D\" or \"-align --all\" with similar path parameters (working directory, batchname, ...)",e);
				}
			}
			try{
				sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + j_align_filename, J_gene , 10 , false , indexed_seqlist , sorted_alignments);
			}
			catch(exception& e){
				return terminate_IGoR_with_error_message("Exception caught while reading J alignments before inference/evaluation. Make sure alignments were carried previously using \"-align --J\" or \"-align --all\" with similar path parameters (working directory, batchname, ...)",e);
			}

			vector<tuple<int,string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

			//create the output directory
			system(&("mkdir " + cl_path +  batchname + "output")[0]);

			if(infer){
				//create inference directory directory
				system(&("mkdir " + cl_path +  batchname + "inference")[0]);
				genmodel.infer_model(sorted_alignments_vec , n_iter_inference , cl_path +  batchname + "inference/" , true , likelihood_thresh_inference , viterbi_inference , proba_threshold_ratio_inference);
			}

			if(evaluate){
				//create evaluate directory
				system(&("mkdir " + cl_path +  batchname + "evaluate")[0]);
				genmodel.infer_model(sorted_alignments_vec , 1 , cl_path +  batchname + "evaluate/" , false , likelihood_thresh_evaluate , viterbi_evaluate , proba_threshold_ratio_evaluate);
			}
		}
		else if(infer and evaluate){
			return terminate_IGoR_with_error_message("Cannot infer and evaluate in a single command, please split in two commands (otherwise the model used to evaluate is ambiguous)");
		} //end infer / evaluate

		if(generate){

			system(&("mkdir " + cl_path +  batchname + "generated")[0]);

			GenModel genmodel(cl_model_parms,cl_model_marginals,cl_counters_list);

			//TODO create generated folder
			string w_err_str;
			if(generate_werr){
				w_err_str = "werr";
			}
			else{
				w_err_str = "noerr";
			}

			//Initialize generated sequences output function
			std::list<std::pair<gen_seq_trans,shared_ptr<void>>> func_data_pairs_list;

			if(gen_output_CDR3_data){
				//Get V and J event
				auto events_map = cl_model_parms.get_events_map();
				shared_ptr<const Rec_Event> v_event_ptr = events_map.at(make_tuple(GeneChoice_t,V_gene,Undefined_side));
				shared_ptr<const Rec_Event> j_event_ptr = events_map.at(make_tuple(GeneChoice_t,J_gene,Undefined_side));

				//Get V and J event positions on the queue
				auto model_queue = cl_model_parms.get_model_queue();
				size_t i=0;
				size_t v_queue_pos;
				size_t j_queue_pos;
				while(not model_queue.empty()){
					if(model_queue.front()->get_name() == v_event_ptr->get_name()){
						v_queue_pos = i;
					}
					else if(model_queue.front()->get_name() == j_event_ptr->get_name()){
						j_queue_pos = i;
					}
					++i;
					model_queue.pop();
				}
				shared_ptr<ostream> outputfile_ptr = shared_ptr<ostream>(new ofstream(cl_path +  batchname + "generated/" +"generated_seqs_" + w_err_str + "_CDR3_info.csv"));
				shared_ptr<void> CDR3_func_data_ptr = shared_ptr<void>(new gen_CDR3_data(v_CDR3_anchors,v_event_ptr->get_realizations_map(),v_queue_pos,
						j_CDR3_anchors,j_event_ptr->get_realizations_map(),j_queue_pos,
						outputfile_ptr));

				func_data_pairs_list.emplace_back(output_CDR3_gen_data,CDR3_func_data_ptr);
			}

			genmodel.generate_sequences(generate_n_seq,generate_werr,
					cl_path +  batchname + "generated/" + gen_filename_prefix +"generated_seqs_" + w_err_str + ".csv",
					cl_path + batchname + "generated/" + gen_filename_prefix +"generated_realizations_" + w_err_str + ".csv",
					func_data_pairs_list,false,gen_random_engine_seed);
		}

	}

	else{
		//Write your custom procedure here

	}



	return EXIT_SUCCESS;

}



