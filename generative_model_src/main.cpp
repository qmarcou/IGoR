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


	//Working directory vars
	bool wd = false;
	string cl_path;

	//Chains vars
	bool chain_provided = false;
	bool has_D = false;

	//Custom genomic templates loading variables
	bool custom_v = false;
	bool custom_d = false;
	bool custom_j = false;

	//Input sequences variables
	bool read_seqs = false;
	bool fasta_seqs = false;
	string input_seqs_file;

	//Genomic templates list and aligns parms
	vector<pair<string,string>> v_genomic;
	vector<pair<string,string>> d_genomic;
	vector<pair<string,string>> j_genomic;

	//Model parms and marginals
	Model_Parms cl_model_parms;
	Model_marginals cl_model_marginals;
	map<size_t,shared_ptr<Counter>> cl_counters_list;

	//Sequence generation parms
	size_t generate_n_seq;
	bool generate_werr = true;
	string generated_batchname = "";

	//Inference parms
	bool viterbi_inference;
	double likelihood_thresh_inference;
	double proba_threshold_ratio_inference;
	size_t n_iter_inference;

	//Sequence evaluation parms
	bool viterbi_evaluate;
	double likelihood_thresh_evaluate;
	double proba_threshold_ratio_evaluate;




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

		else if(string(argv[carg_i]) == "-run_demo"){
			cout<<"running demo code"<<endl;
			run_demo = true;
		}

		else if(string(argv[carg_i]) == "-run_custom"){
			cout<<"running custom code"<<endl;
			run_demo = custom;
		}
		else if(string(argv[carg_i]) == "-align"){
			//Provide a boolean for aligning
			align = true;
		}

		else if( (string(argv[carg_i]) == "-infer") or (string(argv[carg_i]) == "-evaluate")){
			//Provide a boolean for inference
			if(string(argv[carg_i]) == "-infer"){
				infer = true;
			}
			else{
				evaluate = true;
			}

			while(string(argv[carg_i+1]).substr(0,2)=="--"){
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
				}
				else if(string(argv[carg_i]) == "--viterbi"){
					if(string(argv[carg_i]) == "-infer"){
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

					if(string(argv[carg_i]) == "-infer"){
						proba_threshold_ratio_inference = p_ratio;
					}
					else{
						proba_threshold_ratio_evaluate = p_ratio;
					}

				}
				else if(string(argv[carg_i]) == "--N_iter"){
					if(string(argv[carg_i]) == "-infer"){
						try{

						}
						catch(exception& e){
							cout<<"Expected an integer for the number of iterations to perform for the inference, received: \"" + string(argv[carg_i]) + "\""<<endl;
							cout<<"Terminating and throwing exception now..."<<endl;
							throw e;
						}
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
			if(string(argv[carg_i]) == "--alpha"){
				has_D = false;
				v_genomic = read_genomic_fasta("../models/human/tcr_alpha/ref_genome/genomicVs.fasta");
				j_genomic = read_genomic_fasta("../models/human/tcr_alpha/ref_genome/genomicJs.fasta");

			}
			else if(string(argv[carg_i]) == "--beta"){
				has_D = true;
				v_genomic = read_genomic_fasta("../models/human/tcr_beta/ref_genome/genomicVs.fasta");
				d_genomic = read_genomic_fasta("../models/human/tcr_beta/ref_genome/genomicDs.fasta");
				j_genomic = read_genomic_fasta("../models/human/tcr_beta/ref_genome/genomicJs.fasta");
			}
			else if(string(argv[carg_i]) == "--light"){
				cout<<"Support for light chains in command line is not ready yet due to the lack of genomic templates and suitable model"<<endl;
				cout<<"If you wish to use IGoR on light chains please contact us so we can work on incorporating a light chain model to IGoR"<<endl;
				throw invalid_argument("Light chains support does not exist yet for command line");
			}
			else if( (string(argv[carg_i]) == "--heavy_naive") or (string(argv[carg_i]) == "--heavy_memory") ){
				has_D = true;
				v_genomic = read_genomic_fasta("../models/human/bcr_heavy/ref_genome/genomicVs.fasta");
				d_genomic = read_genomic_fasta("../models/human/bcr_heavy/ref_genome/genomicDs.fasta");
				j_genomic = read_genomic_fasta("../models/human/bcr_heavy/ref_genome/genomicJs.fasta");

				if( string(argv[carg_i]) == "--heavy_naive" ){
					//Use a single error rate
				}
				else{
					//Memory

					//TODO infer only \mu for the hypermutation model
				}
			}
			else{
				throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" to specify the chain choice!\n Supported arguments are: --alpha, --beta, --heavy_naive , --heavy_memory , --light");
			}
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

		/*
		 * Sequence generation arguments parsing
		 */
		else if(string(argv[carg_i]) == "-generate"){
			generate = true;
			++carg_i;
			if(string(argv[carg_i]).substr(0,2) == "--"){
				while(string(argv[carg_i]).substr(0,2) == "--"){
					if(string(argv[carg_i]) == "--noerr"){
						generate_werr = false;
					}
					else{
						throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" to specify sequence generation parameters");
					}
				}
				// The number of sequences to be generated has to be given after the arguments
				try{
					generate_n_seq = stoi(string(argv[carg_i]));
				}
				catch(exception& e){
					cout<<"Expected the number of sequences to generate, received: \"" + string(argv[carg_i]) + "\""<<endl;
					cout<<"Terminating and throwing exception now..."<<endl;
					throw e;
				}
			}
			else{
				// or before the other generation arguments
				try{
					generate_n_seq = stoi(string(argv[carg_i]));
				}
				catch(exception& e){
					cout<<"Expected the number of sequences to generate, received: \"" + string(argv[carg_i]) + "\""<<endl;
					cout<<"Terminating and throwing exception now..."<<endl;
					throw e;
				}
				if(string(argv[carg_i]).substr(0,2) == "--"){
					while(string(argv[carg_i]).substr(0,2) == "--"){
						if(string(argv[carg_i]) == "--noerr"){
							generate_werr = false;
						}
						else if(string(argv[carg_i]) == "--name"){
							++carg_i;
							generated_batchname = argv[carg_i];
							++carg_i;
						}
						else{
							throw invalid_argument("Unknown argument \""+string(argv[carg_i])+"\" to specify sequence generation parameters");
						}
					}
				}
			}
		}

		/*
		 * Input sequences argument parsing
		 */
		else if(string(argv[carg_i]) == "-read_seqs"){
			++carg_i;
			string tmp_str = string(argv[carg_i]);
			tmp_str = tmp_str.substr(tmp_str.size() -6 , string::npos );
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

	if(not wd){
		cl_path = "/tmp/";
	}

	cout<<"Working directory set to: \""+cl_path+"\"";




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

			write_indexed_seq_csv(cl_path + "aligns/indexed_sequences.csv",indexed_seqlist);
		}

		if(align){
			vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(cl_path + "aligns/indexed_sequences.csv");

			//Declare substitution matrix used for alignments(nuc44 here)
			double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
			Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

			//Instantiate aligner with substitution matrix declared above and gap penalty of 50

			//Performs V alignments
			cout<<"Performing V alignments...."<<endl;
			Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
			v_aligner.set_genomic_sequences(v_genomic);

			v_aligner.align_seqs(cl_path + "aligns/V_alignments.csv" , indexed_seqlist , 50 , true);

			if(has_D){
				//Performs D alignments if the chain contains a D
				cout<<"Performing D alignments...."<<endl;
				Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
				d_aligner.set_genomic_sequences(d_genomic);

				d_aligner.align_seqs(cl_path + "aligns/D_alignments.csv",indexed_seqlist,0,false);
			}

			cout<<"Performing J alignments...."<<endl;
			Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
			j_aligner.set_genomic_sequences(j_genomic);

			j_aligner.align_seqs(cl_path + "aligns/D_alignments.csv",indexed_seqlist,10,true); //TODO allow for setting offsets bounds
		}

		if(infer xor evaluate){

			GenModel genmodel(cl_model_parms,cl_model_marginals,cl_counters_list);

			//Reading alignments
			vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(cl_path + "aligns/indexed_sequences.csv");

			unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/V_alignments.csv", V_gene , 55 , false , indexed_seqlist  );
			if(has_D){
				sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/D_alignments.csv", D_gene , 35 , false , indexed_seqlist , sorted_alignments);
			}
			sorted_alignments = read_alignments_seq_csv_score_range(cl_path + "aligns/J_alignments.csv", J_gene , 10 , false , indexed_seqlist , sorted_alignments);

			vector<tuple<int,string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

			//create the output directory
			system(&("mkdir " + cl_path + "output")[0]);

			if(infer){
				//create inference directory directory
				system(&("mkdir " + cl_path + "inference")[0]);
				genmodel.infer_model(sorted_alignments_vec , n_iter_inference , cl_path + "inference/" , true , likelihood_thresh_inference , viterbi_inference , proba_threshold_ratio_inference);
			}

			if(evaluate){
				//create evaluate directory
				system(&("mkdir " + cl_path + "evaluate")[0]);
				genmodel.infer_model(sorted_alignments_vec , 1 , cl_path + "evaluate/" , false , likelihood_thresh_evaluate , viterbi_evaluate , proba_threshold_ratio_evaluate);
			}
		}
		else if(infer and evaluate){
			cout<<"Cannot infer and evaluate in a single command, please split in two commands (otherwise the model used to evaluate is ambiguous)"<<endl;
		}

		if(generate){
			GenModel genmodel(cl_model_parms,cl_model_marginals,cl_counters_list);

			//TODO create generated folder
			string w_err_str;
			if(generate_werr){
				w_err_str = "werr";
			}
			else{
				w_err_str = "noerr";
			}

			genmodel.generate_sequences(generate_n_seq,generate_werr,cl_path + "generated/" +  generated_batchname+"generated_seqs_" + w_err_str + ".csv",cl_path + "generated/" +  generated_batchname+"generated_realizations_" + w_err_str + ".csv");
		}

	}

	else{
		//Write your custom procedure here

	}



	return 0;


}



