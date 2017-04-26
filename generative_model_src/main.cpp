/*
 * main.cpp
 *
 *  Created on: Jan 12, 2015
 *      Author: quentin
 */

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
#include "Utils.h"
#include <chrono>
#include<set>


using namespace std;


int main(int argc , char* argv[]){

	//Command line argument iterator
	size_t carg_i = 1;
	//cout<<argv[argc-1]<<endl;

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

	//Input sequences variables
	bool read_seqs = false;
	bool fasta_seqs = false;
	string input_seqs_file;

	//Genomic templates list and aligns parms
	vector<pair<string,string>> v_genomic;
	vector<pair<string,string>> d_genomic;
	vector<pair<string,string>> j_genomic;

	//Model parms and marginals
	bool load_last_inferred_parms = false;
	bool custom_cl_parms = false;
	Model_Parms cl_model_parms;
	Model_marginals cl_model_marginals;
	map<size_t,shared_ptr<Counter>> cl_counters_list;

	//Sequence generation parms
	size_t generate_n_seq;
	bool generate_werr = true;

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
	double heavy_pen_nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
	Matrix<double> heavy_pen_nuc44_sub_matrix(4,4,heavy_pen_nuc44_vect);

		//V alignment vars
		bool align_v = false;
		string v_align_filename = "V_alignments.csv";
		double v_align_thresh_value = 50.0;
		Matrix<double> v_subst_matrix = heavy_pen_nuc44_sub_matrix;
		double v_gap_penalty = 50.0;
		bool v_best_only = true;
		int v_left_offset_bound = INT16_MIN;
		int v_right_offset_bound = INT16_MAX;

		//D alignment vars
		bool align_d = false;
		string d_align_filename = "D_alignments.csv";
		double d_align_thresh_value = 15.0;
		Matrix<double> d_subst_matrix = heavy_pen_nuc44_sub_matrix;
		double d_gap_penalty = 50.0;
		bool d_best_only = false;
		int d_left_offset_bound = INT16_MIN;
		int d_right_offset_bound = INT16_MAX;

		//J alignment vars
		bool align_j = false;
		string j_align_filename = "J_alignments.csv";
		double j_align_thresh_value = 15.0;
		Matrix<double> j_subst_matrix = heavy_pen_nuc44_sub_matrix;
		double j_gap_penalty = 50.0;
		bool j_best_only = true;
		int j_left_offset_bound = INT16_MIN;
		int j_right_offset_bound = INT16_MAX;


	while(carg_i<argc){

		//Command line argument setting the number of threads
		if(string(argv[carg_i]) == string("-threads")){
			omp_set_num_threads(std::stoi(argv[++carg_i]));
			cout<<"Setting number of threads to: "<<argv[carg_i]<<endl;
		}

		//Command line to redirect the standard output to a file
		else if(string(argv[carg_i]) == string("-stdout_f")){
			cout<<"Redirecting output to file: "<<argv[carg_i+1]<<endl;
			freopen(argv[++carg_i],"a+",stdout);
		}

		else if(string(argv[carg_i]) == string("-batch")){
			++carg_i;
			batchname = argv[carg_i];
			if(batchname[batchname.size()-1] != '_'){
				batchname.append("_");
			}
			cout<<"Batch name set to: "<<batchname<<endl;
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
					throw invalid_argument("Unknown gene argument \"" + string(argv[carg_i]) +"\" to set genomic templates");
				}
			}

			if((not custom_v)
					and (not custom_d)
					and (not custom_j)){
				throw runtime_error("No gene argument was passed after -set_genomic");
			}
		}

		else if(string(argv[carg_i]) == "-set_custom_model"){
			custom_cl_parms = true;
			++carg_i;
			cl_model_parms.read_model_parms(string(argv[carg_i]));
			++carg_i;
			cl_model_marginals = Model_marginals(cl_model_parms);
			cl_model_marginals.txt2marginals(string(argv[carg_i]),cl_model_parms);

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
			cout<<"running demo code"<<endl;
			run_demo = true;
		}

		else if(string(argv[carg_i]) == "-run_custom"){
			cout<<"running custom code"<<endl;
			custom = true;
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

				bool best_only_provided = false;
				bool best_only;

				bool offset_bounds_provided = false;
				int left_offset_bound;
				int right_offset_bound;

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
								cout<<"Expected a float for the alignment score threshold, received: \"" + string(argv[carg_i]) + "\""<<endl;
								cout<<"Terminating and throwing exception now..."<<endl;
								throw e;
							}
							thresh_provided = true;
						}
						else if(string(argv[carg_i]) == "---matrix"){
							//Read the substitution matrix from a file
							//TODO code this
							throw invalid_argument("Support for a changing the substitution matrix from a file is not coded yet");
						}
						else if(string(argv[carg_i]) == "---gap_penalty"){
							//Read the gap penalty
							++carg_i;
							try{
								gap_penalty = stod(string(argv[carg_i]));
							}
							catch (exception& e) {
								cout<<"Expected a float for the alignment gap penalty, received: \"" + string(argv[carg_i]) + "\""<<endl;
								cout<<"Terminating and throwing exception now..."<<endl;
								throw e;
							}
							gap_penalty_provided = true;
						}
						else if(string(argv[carg_i]) == "---best_only"){
							//Set the best alignment only boolean
							++carg_i;
							if(string(argv[carg_i]) == "true"){
								best_only = true;
							}
							else if(string(argv[carg_i]) == "false"){
								best_only = false;
							}
							else{
								throw invalid_argument("Unkown argument received\"" + string(argv[carg_i]) + "\"to set best alignment only boolean, existing values are: true, false");
							}
							best_only_provided = true;

						}
						else if(string(argv[carg_i]) == "---offset_bounds"){
							//Read the offset bounds
							++carg_i;
							try{
								left_offset_bound = stoi(string(argv[carg_i]));
							}
							catch (exception& e) {
								cout<<"Expected an integer for the left offset bound, received: \"" + string(argv[carg_i]) + "\""<<endl;
								cout<<"Terminating and throwing exception now..."<<endl;
								throw e;
							}

							++carg_i;
							try{
								right_offset_bound = stoi(string(argv[carg_i]));
							}
							catch (exception& e) {
								cout<<"Expected an integer for the right offset bound, received: \"" + string(argv[carg_i]) + "\""<<endl;
								cout<<"Terminating and throwing exception now..."<<endl;
								throw e;
							}
							offset_bounds_provided = true;
						}
						else{
							throw invalid_argument("Unknown parameter\"" + string(argv[carg_i]) +"\" for gene " + gene_str_val + " in -align " );
						}
					}
				}
				else{
					throw invalid_argument("Unknown gene specification\"" + string(argv[carg_i]) + "\"for -align");
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

					if(best_only_provided){
						v_best_only = best_only;
					}

					if(offset_bounds_provided){
						v_left_offset_bound = left_offset_bound;
						v_right_offset_bound = right_offset_bound;
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

					if(best_only_provided){
						d_best_only = best_only;
					}

					if(offset_bounds_provided){
						d_left_offset_bound = left_offset_bound;
						d_right_offset_bound = right_offset_bound;
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

					if(best_only_provided){
						j_best_only = best_only;
					}

					if(offset_bounds_provided){
						j_left_offset_bound = left_offset_bound;
						j_right_offset_bound = right_offset_bound;
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
						cout<<"Expected a float for the likelihood threshold, received: \"" + string(argv[carg_i]) + "\""<<endl;
						cout<<"Terminating and throwing exception now..."<<endl;
						throw e;
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
						cout<<"Expected a float for the probability ratio threshold, received: \"" + string(argv[carg_i]) + "\""<<endl;
						cout<<"Terminating and throwing exception now..."<<endl;
						throw e;
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
							cout<<"Expected an integer for the number of iterations to perform for the inference, received: \"" + string(argv[carg_i]) + "\""<<endl;
							cout<<"Terminating and throwing exception now..."<<endl;
							throw e;
						}
					}
					else{
						throw invalid_argument("Invalid argument \"--N_iter\" for -evaluate");
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
				else if(string(argv[carg_i]) == "--subsample"){
					subsample_seqs = true;
					++carg_i;
					try{
						n_subsample_seqs = stoi(argv[carg_i]);
					}
					catch(exception& e){
						cout<<"Expected an integer for the number of sequences to subsample, received: \"" + string(argv[carg_i]) + "\""<<endl;
						cout<<"Terminating and throwing exception now..."<<endl;
						throw e;
					}
				}
				else{
					throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" to specify inference/evaluate parameters");
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
				cout<<"Chain parameter set to: "<<chain_arg_str<<endl;
			}
			else{
				throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" to specify the chain choice!\n Supported arguments are: alpha, beta, heavy_naive , heavy_memory , light");
			}
		}

		else if(string(argv[carg_i]) == "-species"){
			//TODO add a check on the existence of the species
			species_provided = true;
			++carg_i;
			species_str = string(argv[carg_i]);
			cout<<"Species parameter set to: "<<species_str<<endl;
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
				throw invalid_argument("Working directory needs to be set before declaring the counters, please re-order the arguments");
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
						cout<<"Expected the number of scenarios to be recorded by the best scenario counter, received: \"" + string(argv[carg_i]) + "\""<<endl;
						cout<<"Terminating and throwing exception now..."<<endl;
						throw e;
					}

					if(n_record_scenarios<=0){
						throw invalid_argument("Number of scenarios to be recorded must be greater than zero");
					}

					shared_ptr<Counter>best_sc_ptr(new Best_scenarios_counter(10 , cl_path + "output/" ,true));
					cl_counters_list.emplace(cl_counters_list.size(),best_sc_ptr);
				}
				else if(string(argv[carg_i]) == "--coverage"){
					Gene_class chosen_gc;
					++carg_i;
					try{
						chosen_gc = str2GeneClass(string(argv[carg_i]));
					}
					catch(exception& e){
						throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" to specify coverage target!\n Supported arguments are: V_gene, VD_genes, D_gene, DJ_gene, VJ_gene, J_gene, VDJ_genes");
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
				cout<<"Expected the number of sequences to generate, received: \"" + string(argv[carg_i]) + "\""<<endl;
				cout<<"Terminating and throwing exception now..."<<endl;
				throw e;
			}

			while( (carg_i+1<argc)
					and (string(argv[carg_i+1]).size()>2)
					and string(argv[carg_i+1]).substr(0,2) == "--"){

				++carg_i;

				if(string(argv[carg_i]) == "--noerr"){
					generate_werr = false;
				}
				else{
					throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" to specify sequence generation parameters");
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
			string tmp_str = input_seqs_file.substr(input_seqs_file.size() -6 , string::npos );
			transform(tmp_str.begin(),tmp_str.end(),tmp_str.begin(),::tolower);
			if( tmp_str == ".fasta" ){
				fasta_seqs = true;
				cout<<"FASTA extension detected for the input sequence file"<<endl;
			}
			else{
				cout<<"No FASTA extension detected for the input sequence file assuming a text file without header"<<endl;
			}
		}

		//If the argument does not correspond to any previous section throw an exception
		else{
			throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" ");
		}

		//Read the next command line argument
		++carg_i;
	}

	//Make sure the working directory is set somewhere before performing any action
	if(not wd){
		cl_path = "/tmp/";
	}
	cout<<"Working directory set to: \""+cl_path+"\""<<endl;

	//Check that both species and chain have been provided
	if(chain_provided xor species_provided){
		throw ("Both species and chain must be provided when using a predefined model!");
	}


	if(chain_provided){
		if(chain_arg_str == "alpha"){
			has_D = false;
			chain_path_str = "tcr_alpha";
			v_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			j_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");

		}
		else if(chain_arg_str == "beta"){
			has_D = true;
			chain_path_str = "tcr_beta";
			v_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			d_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicDs.fasta");
			j_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");
		}
		else if(chain_arg_str == "light"){
			cout<<"Support for light chains in command line is not ready yet due to the lack of genomic templates and suitable model"<<endl;
			cout<<"If you wish to use IGoR on light chains please contact us so we can work on incorporating a light chain model to IGoR"<<endl;
			throw invalid_argument("Light chains support does not exist yet for command line");
		}
		else if( (chain_arg_str == "heavy_naive") or (chain_arg_str == "heavy_memory") ){
			has_D = true;
			chain_path_str = "bcr_heavy";
			v_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicVs.fasta");
			d_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicDs.fasta");
			j_genomic = read_genomic_fasta("../models/"+species_str+"/"+chain_path_str+"/ref_genome/genomicJs.fasta");

			if( chain_arg_str == "heavy_naive" ){
				//Use a single error rate
			}
			else{
				//Memory

				//TODO infer only \mu for the hypermutation model
			}
		}
	}

	//Read custom genomic templates if some custom ones were specified
	if(custom_v){
		v_genomic = read_genomic_fasta(custom_v_path);
	}
	if(custom_d){
		has_D = true;
		d_genomic = read_genomic_fasta(custom_d_path);
	}
	if(custom_j){
		j_genomic = read_genomic_fasta(custom_j_path);
	}

	//Make sure passed arguments are unambiguous
	if(custom_cl_parms and load_last_inferred_parms){
		throw invalid_argument("Setting a custom model and loading the last inferred model in the same command is ambiguous!");
	}

	//Load last inferred model
	if(load_last_inferred_parms){
		cout<<"Loading last inferred model..."<<endl;
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
			cout<<"Failed to load last inferred model, please check that the model exists"<<endl;
			cout<<"Throwing exception now..."<<endl;
			throw e;
		}
	}

	/*
	 * Read supplied model parms and marginals
	 * If some custom genomic templates were supplied, we replace the genomic templates contained in the model by the supplied ones
	 * and re-initialize the marginals
	 */
	if( ((not custom_cl_parms) and (not load_last_inferred_parms))
			and (infer or evaluate or generate)){
		cout<<"read some model parms"<<endl;
		cl_model_parms.read_model_parms("../models/"+species_str+"/"+chain_path_str+"/models/model_parms.txt");


		bool any_custom_gene = false;
		unordered_map<tuple<Event_type,Gene_class,Seq_side>,shared_ptr<Rec_Event>> tmp_events_map = cl_model_parms.get_events_map();
		if(custom_v){
			any_custom_gene = true;
			shared_ptr<Rec_Event> v_choice = tmp_events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side));
			shared_ptr<Gene_choice> v_choice_gc = dynamic_pointer_cast<Gene_choice>(v_choice);
			Rec_Event_name former_name = v_choice_gc->get_name();
			v_choice_gc->set_genomic_templates(v_genomic);
			cl_model_parms.update_edge_event_name(former_name,v_choice_gc->get_name());
		}
		if(has_D and custom_d){
			any_custom_gene = true;
			shared_ptr<Rec_Event> d_choice = tmp_events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side));
			shared_ptr<Gene_choice> d_choice_gc = dynamic_pointer_cast<Gene_choice>(d_choice);
			Rec_Event_name former_name = d_choice_gc->get_name();
			d_choice_gc->set_genomic_templates(d_genomic);
			cl_model_parms.update_edge_event_name(former_name,d_choice_gc->get_name());
		}
		if(custom_j){
			any_custom_gene = true;
			shared_ptr<Rec_Event> j_choice = tmp_events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side));
			shared_ptr<Gene_choice> j_choice_gc = dynamic_pointer_cast<Gene_choice>(j_choice);
			Rec_Event_name former_name = j_choice_gc->get_name();
			j_choice_gc->set_genomic_templates(j_genomic);
			cl_model_parms.update_edge_event_name(former_name,j_choice_gc->get_name());
		}

		cl_model_marginals = Model_marginals(cl_model_parms);
		if(any_custom_gene){
			/*
			 * If some custom genomic templates were provided we have replaced the genomic templates contained in the model
			 * Thus other components (e.g deletion profiles) might not match anymore and re inferring a model is necessary
			 * Marginals are initialized with a uniform distribution
			 */

			cl_model_marginals.uniform_initialize(cl_model_parms);
		}
		else{
			cl_model_marginals.txt2marginals("../models/"+species_str+"/"+chain_path_str+"/models/model_marginals.txt",cl_model_parms);
		}
	}

	/*
	 * Once model parms have been read fix the events requested
	 */
	if(infer_only and no_infer){
		throw invalid_argument("Cannot use both \"--infer_only\" and \"--not_infer\" since they are somewhat redundant");
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

		cout<<"Reading genomic templates"<<endl;

		vector<pair<string,string>> v_genomic = read_genomic_fasta( string("../demo/genomicVs_with_primers.fasta"));

		vector<pair<string,string>> d_genomic = read_genomic_fasta( string("../demo/genomicDs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta( string("../demo/genomicJs_all_curated.fasta"));

		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

		//Instantiate aligner with substitution matrix declared above and gap penalty of 50
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
		d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		cout<<"Reading sequences and aligning"<<endl;
		typedef std::chrono::system_clock myclock;
		myclock::time_point begin_time, end_time;

		begin_time = myclock::now();


		vector<pair<const int, const string>> indexed_seqlist = read_txt( string("../demo/murugan_naive1_noncoding_demo_seqs.txt") ); //Could also read a FASTA file <code>read_fasta()<\code> or indexed sequences <code>read_indexed_seq_csv()<\code>


		v_aligner.align_seqs( string("../demo/murugan_naive1_noncoding_demo_seqs") + string("_alignments_V.csv"),indexed_seqlist,50,true,INT16_MIN,-155);
		//v_aligner.write_alignments_seq_csv(path + string("alignments_V.csv") , v_alignments);

		d_aligner.align_seqs(string("../demo/murugan_naive1_noncoding_demo_seqs") + string("_alignments_D.csv"),indexed_seqlist,0,false);
		//d_aligner.write_alignments_seq_csv(path + string("alignments_D.csv") , d_alignments);

		j_aligner.align_seqs(string("../demo/murugan_naive1_noncoding_demo_seqs") + string("_alignments_J.csv"),indexed_seqlist,10,true,42,48);

		end_time= myclock::now();
		chrono::duration<double> elapsed = end_time - begin_time;
		cout<<"Alignments procedure lasted: "<<elapsed.count()<<" seconds"<<endl;
		cout<<"for "<<indexed_seqlist.size()<<" TCRb sequences of 60bp(from murugan and al), against ";
		cout<<v_genomic.size()<<" Vs,"<<d_genomic.size()<<" Ds, and "<<j_genomic.size()<<" Js full sequences"<<endl;

		//unordered_map<int,forward_list<Alignment_data>> j_alignments = j_aligner.align_seqs(indexed_seqlist,10,true,42,48);
		//j_aligner.write_alignments_seq_csv(path + string("alignments_J.csv") , j_alignments);


		write_indexed_seq_csv(string("../demo/murugan_naive1_noncoding_demo_seqs") + string("_indexed_seq.csv") , indexed_seqlist);



		cout<<"Construct the model"<<endl;
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

		cout<<"Write and read back the model"<<endl;
		//Write the model_parms into a file
		parms.write_model_parms(string("../demo/demo_write_model_parms.txt"));

		//Write the marginals into a file
		model_marginals.write2txt(string("../demo/demo_write_model_marginals.txt"),parms);

		//Read a model and marginal pair
		Model_Parms read_model_parms;
		read_model_parms.read_model_parms(string("../demo/demo_write_model_parms.txt"));
		Model_marginals read_model_marginals(read_model_parms);
		read_model_marginals.txt2marginals(string("../demo/demo_write_model_marginals.txt"),read_model_parms);

		//Instantiate a Counter
		map<size_t,shared_ptr<Counter>> counters_list;
		 //Collect gene coverage and errors
		shared_ptr<Counter> coverage_counter_ptr(new Coverage_err_counter("../demo/run_demo/",VJ_genes,1,false,false));
		counters_list.emplace(0,coverage_counter_ptr);

		 //Collect 10 best scenarios per sequence during the last iteration
		shared_ptr<Counter>best_sc_ptr(new Best_scenarios_counter(10 , "../demo/run_demo/" ,true ));
		counters_list.emplace(1,best_sc_ptr);

		 //Collect sequence generation probability during last iteration
		shared_ptr<Counter> pgen_counter_ptr(new Pgen_counter ("../demo/run_demo/"));
		counters_list.emplace(2,pgen_counter_ptr);

		//Instantiate the high level GenModel class
		//This class allows to make most useful high level operations(model inference/Pgen computation , sequence generation)
		GenModel gen_model(read_model_parms,read_model_marginals,counters_list);


		//Inferring a model

		//Read alignments
		//vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path+ string(argv[2]) + string("indexed_seq.csv"));
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(string("../demo/murugan_naive1_noncoding_demo_seqs") + string("_alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
		sorted_alignments = read_alignments_seq_csv_score_range(string("../demo/murugan_naive1_noncoding_demo_seqs") + string("_alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15
		sorted_alignments = read_alignments_seq_csv_score_range(string("../demo/murugan_naive1_noncoding_demo_seqs") + string("_alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<tuple<int,string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

		//Infer the model
		cout<<"Infer model"<<endl;

		begin_time = myclock::now();

		gen_model.infer_model(sorted_alignments_vec , 20 , string("../demo/run_demo/") , true ,1e-35,0.0001);

		end_time= myclock::now();
		elapsed = end_time - begin_time;
		cout<<"Model inference procedure lasted: "<<elapsed.count()<<" seconds"<<endl;


		//Generate sequences
		cout<<"Generate sequences"<<endl;
		auto generated_seq =  gen_model.generate_sequences(100,false);//Without errors
		gen_model.write_seq_real2txt(string("../demo/generated_seqs_indexed_demo.csv"), string("../demo/generated_seqs_realizations_demo.csv") , generated_seq);//Member function will be changed

	}

	else if (not custom){
		//Execute code dictated by command line arguments
		if(read_seqs){
			vector<pair<const int, const string>> indexed_seqlist;
			if(fasta_seqs){
				indexed_seqlist = read_fasta(input_seqs_file);
			}
			else{
				indexed_seqlist = read_txt(input_seqs_file);
			}
			//create the directory
			system(&("mkdir " + cl_path + "aligns")[0]);

			write_indexed_seq_csv(cl_path + "aligns/" + batchname + "indexed_sequences.csv",indexed_seqlist);
		}

		if(align){
			vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(cl_path + "aligns/" + batchname + "indexed_sequences.csv");

			if(align_v){
				//Performs V alignments
				cout<<"Performing V alignments...."<<endl;
				Aligner v_aligner = Aligner(v_subst_matrix , v_gap_penalty , V_gene);
				v_aligner.set_genomic_sequences(v_genomic);

				v_aligner.align_seqs(cl_path + "aligns/" +  batchname + v_align_filename , indexed_seqlist , v_align_thresh_value , v_best_only , v_left_offset_bound , v_right_offset_bound);
			}


			if(has_D and align_d){
				//Performs D alignments if the chain contains a D
				cout<<"Performing D alignments...."<<endl;
				Aligner d_aligner = Aligner(d_subst_matrix , d_gap_penalty , D_gene);
				d_aligner.set_genomic_sequences(d_genomic);

				d_aligner.align_seqs(cl_path + "aligns/" +  batchname + d_align_filename ,indexed_seqlist, d_align_thresh_value , d_best_only , d_left_offset_bound , d_right_offset_bound);
			}

			if(align_j){
				cout<<"Performing J alignments...."<<endl;
				Aligner j_aligner (j_subst_matrix , j_gap_penalty , J_gene);
				j_aligner.set_genomic_sequences(j_genomic);

				j_aligner.align_seqs(cl_path + "aligns/" +  batchname + j_align_filename , indexed_seqlist, j_align_thresh_value , j_best_only , j_left_offset_bound , j_right_offset_bound);
			}

		}

		if(infer xor evaluate){

			GenModel genmodel(cl_model_parms,cl_model_marginals,cl_counters_list);

			//Reading alignments
			vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(cl_path + "aligns/" + batchname + "indexed_sequences.csv");

			if(subsample_seqs){
				indexed_seqlist = sample_indexed_seq(indexed_seqlist,n_subsample_seqs);
			}

			unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + v_align_filename, V_gene , 55 , false , indexed_seqlist  );
			if(has_D){
				sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + d_align_filename, D_gene , 35 , false , indexed_seqlist , sorted_alignments);
			}
			sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/" +  batchname + j_align_filename, J_gene , 10 , false , indexed_seqlist , sorted_alignments);

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
			cout<<"Cannot infer and evaluate in a single command, please split in two commands (otherwise the model used to evaluate is ambiguous)"<<endl;
		}

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

			genmodel.generate_sequences(generate_n_seq,generate_werr,cl_path +  batchname + "generated/" +"generated_seqs_" + w_err_str + ".csv",cl_path + "generated/" + batchname +"generated_realizations_" + w_err_str + ".csv");
		}

	}

	else{
		//Write your custom procedure here

	}



	return 0;


}



