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
#include <chrono>


using namespace std;


int main(int argc , char* argv[]){

	size_t carg_i = 1;
	cout<<argv[argc-1]<<endl;
	bool run_demo = false;
	bool align = false;
	bool infer = false;

	while(carg_i<argc){

		//Command line argument setting the number of threads
		if(string(argv[carg_i]) == string("-threads")){
			omp_set_num_threads(std::stoi(argv[++carg_i]));
			cout<<"Setting number of threads to: "<<argv[carg_i]<<endl;
		}

		//Command line to redirect the standard output to a file
		if(string(argv[carg_i]) == string("-stdout_f")){
			cout<<"Redirecting output to file: "<<argv[carg_i+1]<<endl;
			freopen(argv[++carg_i],"a+",stdout);
		}

		if(string(argv[carg_i]) == "-run_demo"){
			cout<<"running demo code"<<endl;
			run_demo = true;
		}

		if(string(argv[carg_i]) == "-align"){
			//Provide a boolean for aligning
			align = true;
		}

		if(string(argv[carg_i]) == "-infer"){
			//Provide a boolean for infering
			infer = true;
		}

		//Read the next command line argument
		++carg_i;


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

		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

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


	else{
		//Write your custom procedure here

		string path = "/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/run_no_d/shm_test/";



		//Read the alpha chain model of BCRs
		Model_Parms alpha_model_parms;
		alpha_model_parms.read_model_parms("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/run_no_d/iteration_20_parms.txt");
		Model_marginals alpha_model_marginals(alpha_model_parms);
		alpha_model_marginals.txt2marginals("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/run_no_d/iteration_20.txt",alpha_model_parms);


		//Instantiate a Hypermutation model and set it as the new error rate
		Hypermutation_global_errorrate shm_err_rate (3,V_gene,VDJ_genes,.01);
		//shm_err_rate.generate_random_contributions(10);
		alpha_model_parms.set_error_ratep(&shm_err_rate);
		//alpha_model_parms.write_model_parms("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/run_no_d/shm_test/test.txt");


/*		//Read back the model
		Model_Parms alpha_read_parms;
		alpha_read_parms.read_model_parms("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/run_no_d/shm_test/test.txt");

		//introduce errors in an out of frame naive sequences sample
		auto naive_indexed_seq = read_indexed_csv("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/AJ_Naive_noncoding_indexed_seq.csv");

		//Create seed for random generator
		//create a seed from timer
		typedef std::chrono::high_resolution_clock myclock;
		myclock::time_point time = myclock::now();
		myclock::duration dur = myclock::time_point::max() - time;

		unsigned time_seed = dur.count();
		//Instantiate random number generator
		default_random_engine generator =  default_random_engine(time_seed);

		vector<pair<const int,const string>> mutated_indexed_seq;
		for(auto iter = naive_indexed_seq.begin() ; iter != naive_indexed_seq.end() ; ++iter){
			string str_copy = (*iter).second;
			shm_err_rate.generate_errors(str_copy,generator);
			mutated_indexed_seq.emplace_back((*iter).first,str_copy);
		}

		write_indexed_seq_csv(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/run_no_d/shm_test/mutated_seqs_test.csv"),mutated_indexed_seq);*/

		//infer back a model keeping everything fixed but the hypermutation error rate
		alpha_model_parms.set_fixed_all_events(true);

/*		//align the sequences
        vector<pair<string,string>> v_genomic = read_genomic_fasta(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicVs.fasta"));

        vector<pair<string,string>> j_genomic = read_genomic_fasta(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicJs.fasta"));



		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-4,-4,-4,1,-4,-4,1,-4,1,-4,-1,-1,-1,-2,-4,5,-4,-4,-4,1,-4,1,1,-4,-1,-4,-1,-1,-2,-4,-4,5,-4,1,-4,1,-4,1,-4,-1,-1,-4,-1,-2,-4,-4,-4,5,-4,1,1,-4,-4,1,-1,-1,-1,-4,-2,1,-4,1,-4,-1,-4,-2,-2,-2,-2,-3,-1,-3,-1,-1,-4,1,-4,1,-4,-1,-2,-2,-2,-2,-1,-3,-1,-3,-1,-4,-4,1,1,-2,-2,-1,-4,-2,-2,-1,-1,-3,-3,-1,1,1,-4,-4,-2,-2,-4,-1,-2,-2,-3,-3,-1,-1,-1,-4,1,1,-4,-2,-2,-2,-2,-1,-4,-1,-3,-3,-1,-1,1,-4,-4,1,-2,-2,-2,-2,-4,-1,-3,-1,-1,-3,-1,-4,-1,-1,-1,-3,-1,-1,-3,-1,-3,-1,-2,-2,-2,-1,-1,-4,-1,-1,-1,-3,-1,-3,-3,-1,-2,-1,-2,-2,-1,-1,-1,-4,-1,-3,-1,-3,-1,-3,-1,-2,-2,-1,-2,-1,-1,-1,-1,-4,-1,-3,-3,-1,-1,-3,-2,-2,-2,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		Matrix<double> nuc44_sub_matrix(15,15,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty

		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);


		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		auto mutated_indexed_seq_sample = sample_indexed_seq(mutated_indexed_seq,40000);

		v_aligner.align_seqs( path + "V_alignments_non_coding_mutated.csv" , mutated_indexed_seq_sample,50,true,INT16_MIN,-145);

		j_aligner.align_seqs(path + "J_alignments_non_coding_mutated.csv" , mutated_indexed_seq_sample,10,true,89,104);*/

		//read alignments
		auto mutated_seqs_sample = read_indexed_csv(path + "mutated_indexed_seq_40000_sample_sub.csv");

		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(path + "V_alignments_non_coding_mutated.csv", V_gene , 55 , false , mutated_seqs_sample  );//40//35

		sorted_alignments = read_alignments_seq_csv_score_range(path + "J_alignments_non_coding_mutated.csv", J_gene , 10 , false , mutated_seqs_sample , sorted_alignments);//30//20

		//infer back the model
		GenModel genmodel(alpha_model_parms,alpha_model_marginals);


		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

		genmodel.infer_model(sorted_alignments_vec , 10 , path+"/run_test/" ,1e-100,0.001);

	}



	return 0;


}



