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


using namespace std;


int main(int argc , char* argv[]){

	omp_set_num_threads(std::stoi(argv[5]));

	//I/O stream to ask the directory for sequences,output etc
		//I/O stream to ask the output file name
		//I/O stream to ask number of iterations/likelihood to reach
		//I/O stream to load existing model?

	if(false){
	/*


		 // Deal with anand data

/*		string path;
		cout<<"enter working directory path:"<<endl;
		cin>>path;*/
		//string path("/home/marcou/tcr_murugan/inference_data/");
		//string path("/home/marcou/tcr_murugan/gen_seq_inf/");
		//string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/tcr_murugan/inference_data/");
		//string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/tcr_murugan/gen_seq_inf/");
		//string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/tcr_murugan/inference_data_subs/");
		//string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/tcr_murugan/one_seq_comp/");
		//string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/tcr_murugan/non_zero_ins_seqs/");
		//string path("/users/marcou/PhD/tcr_murugan/inference_data/inference_data/");
		string path (argv[1]);
		if (path[path.size()-1] != '/'){
			path+="/";
		}



		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("genomicVs_with_primers.fasta"));

		vector<pair<string,string>> d_genomic = read_genomic_fasta(path + string("genomicDs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("genomicJs_all_curated.fasta"));

/*		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner d_aligner = Aligner(nuc44_sub_matrix , 50);
		d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50);
		j_aligner.set_genomic_sequences(j_genomic);

		vector<pair<const int, const string>> indexed_seqlist = read_txt(path + string(argv[3]));

		unordered_map<int,forward_list<Alignment_data>> v_alignments = v_aligner.align_seqs(indexed_seqlist,50,true,INT16_MIN,-155);
		v_aligner.write_alignments_seq_csv(path + string("alignments_V.csv") , v_alignments);

		unordered_map<int,forward_list<Alignment_data>> d_alignments = d_aligner.align_seqs(indexed_seqlist,0,false);
		d_aligner.write_alignments_seq_csv(path + string("alignments_D.csv") , d_alignments);

		unordered_map<int,forward_list<Alignment_data>> j_alignments = j_aligner.align_seqs(indexed_seqlist,10,true,42,48);
		j_aligner.write_alignments_seq_csv(path + string("alignments_J.csv") , j_alignments);
		j_aligner.write_indexed_seq_csv(path + string("indexed_seq.csv") , indexed_seqlist);*/

		//Construct tcr model
		Gene_choice v_choice(V_gene,v_genomic);
		v_choice.set_nickname("v_choice");
		v_choice.set_priority(7);
		Gene_choice d_choice(D_gene,d_genomic);
		d_choice.set_nickname("d_gene");
		d_choice.set_priority(8);
		Gene_choice j_choice(J_gene,j_genomic);
		j_choice.set_nickname("j_choice");
		j_choice.set_priority(6);

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
		parms.add_edge(&d_choice,&j_choice);

		Model_marginals model_marginals(parms);
		model_marginals.uniform_initialize(parms);

		Single_error_rate error_rate(stod(string(argv[4])));

		parms.set_error_ratep(&error_rate);

		GenModel gen_model(parms,model_marginals);

		//Read alignments
		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("indexed_seq.csv"));
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
		sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15
		sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

		cout<<"num_procs: "<<omp_get_num_procs()<<endl;
		cout<<"num_threads: "<<omp_get_num_threads()<<endl;


		gen_model.infer_model(sorted_alignments_vec , 10 , path+string("run") + string(argv[2]) + string("/"));

	}

	if(false){
		//Align Harlan BCR sequences
		vector<pair<string,string>> v_genomic = read_genomic_fasta(string("/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicVs.fasta"));

		vector<pair<string,string>> d_genomic = read_genomic_fasta(string("/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicDs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(string("/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicJs.fasta"));

		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-4,-4,-4,1,-4,-4,1,-4,1,-4,-1,-1,-1,-2,-4,5,-4,-4,-4,1,-4,1,1,-4,-1,-4,-1,-1,-2,-4,-4,5,-4,1,-4,1,-4,1,-4,-1,-1,-4,-1,-2,-4,-4,-4,5,-4,1,1,-4,-4,1,-1,-1,-1,-4,-2,1,-4,1,-4,-1,-4,-2,-2,-2,-2,-3,-1,-3,-1,-1,-4,1,-4,1,-4,-1,-2,-2,-2,-2,-1,-3,-1,-3,-1,-4,-4,1,1,-2,-2,-1,-4,-2,-2,-1,-1,-3,-3,-1,1,1,-4,-4,-2,-2,-4,-1,-2,-2,-3,-3,-1,-1,-1,-4,1,1,-4,-2,-2,-2,-2,-1,-4,-1,-3,-3,-1,-1,1,-4,-4,1,-2,-2,-2,-2,-4,-1,-3,-1,-1,-3,-1,-4,-1,-1,-1,-3,-1,-1,-3,-1,-3,-1,-2,-2,-2,-1,-1,-4,-1,-1,-1,-3,-1,-3,-3,-1,-2,-1,-2,-2,-1,-1,-1,-4,-1,-3,-1,-3,-1,-3,-1,-2,-2,-1,-2,-1,-1,-1,-1,-4,-1,-3,-3,-1,-1,-3,-2,-2,-2,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		Matrix<double> nuc44_sub_matrix(15,15,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty

		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

	 //   Aligner d_aligner = Aligner(nuc44_sub_matrix , 50);
	   // d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		vector<pair<const int, const string>> indexed_seqlist = read_txt("/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/MEMORY_04-AJ-M_001-025/concatenated_04-AJ-M_unique_noncoding_seqonly.txt");

		string path = "/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/MEMORY_04-AJ-M_001-025/";

		v_aligner.align_seqs( path + "V_alignments_non_coding.csv" , indexed_seqlist,50,true,INT16_MIN,-145);
		//v_aligner.write_alignments_seq_csv(path + string("alignments_V.csv") , v_alignments);

		// unordered_map<int,forward_list<Alignment_data>> d_alignments = d_aligner.align_seqs(indexed_seqlist,0,false);
		// d_aligner.write_alignments_seq_csv(path + string("alignments_D.csv") , d_alignments);

		j_aligner.align_seqs(path + "J_alignments_non_coding.csv" , indexed_seqlist,10,true,89,104);
		//j_aligner.write_alignments_seq_csv(path + string("alignments_J.csv") , j_alignments);
		write_indexed_seq_csv(path + string("indexed_seq_non_coding.csv") , indexed_seqlist);
	}

	if(false){

		//string batch("TV_F");
		string batch(argv[1]);
		cout<<"processing "<<batch<<endl;
/*		string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/old_twins_alpha/");

		Model_Parms parms;
		parms.read_model_parms(path + batch + "/run_lowprobthresh/iteration_20_parms.txt");

		Model_marginals marginals(parms);
		marginals.txt2marginals(path + batch + "/run_lowprobthresh/iteration_20.txt",parms);*/

		string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/new_twins_beta/");

		Model_Parms parms;
		parms.read_model_parms(path + batch + "/run_clust/iteration_20_parms.txt");

		Model_marginals marginals(parms);
		marginals.txt2marginals(path + batch + "/run_clust/iteration_20.txt",parms);


		GenModel gen_model(parms,marginals);


		gen_model.generate_sequences(4000000,false,path + batch + "/gen_seq_noerr_Yzh/" + batch + "_for_Yzh2_4000000_genseqs_noerr.csv",path + batch + "/gen_seq_noerr_Yzh/" + batch + "_for_Yzh2_4000000_genseqs_noerr_realizations.csv");

		/*//Cut the sequences
		for(auto seq = generated_seq.begin() ; seq != generated_seq.end() ; seq++){
			string tmp = (*seq).first.substr((*seq).first.size()-101,100);
			(*seq).first = tmp;
		}*/

		marginals.write2txt(path + batch + "/gen_seq_noerr_Yzh/" + batch + "_generation_marginals2.txt",parms);
		parms.write_model_parms(path + batch + "/gen_seq_noerr_Yzh/" + batch + "_generation_model2.txt");

		//gen_model.write_seq_real2txt(path + batch + "/gen_seq_werr/" + batch + "_1000000_genseqs_werr.csv",path + batch + "/gen_seq_werr/" + batch + "_1000000_genseqs_werr_realizations.csv",generated_seq);

	}

	if(false){
		//Compute out of frame probas for Misha twin generated

		string path (argv[1]);
		if (path[path.size()-1] != '/'){
			path+="/";
		}

		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("../genomicVs_with_primers.fasta"));

		vector<pair<string,string>> d_genomic = read_genomic_fasta(path + string("../genomicDs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("../genomicJs_all_curated.fasta"));

		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
		d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		vector<pair<const int, const string>> indexed_seqlist = read_txt(path + string(argv[3]));

		v_aligner.align_seqs(path + string("non_prod_alignments_V.csv"),indexed_seqlist,50,true,INT16_MIN,-115);


		d_aligner.align_seqs(path + string("non_prod_alignments_D.csv"),indexed_seqlist,0,false);


		j_aligner.align_seqs(path + string("non_prod_alignments_J.csv"),indexed_seqlist,10,true,45,57);
		write_indexed_seq_csv(path + string("indexed_seq.csv") , indexed_seqlist);

		Model_Parms parms;
		parms.read_model_parms(path + "generation_model.txt");
		Model_marginals marginals(parms);
		marginals.txt2marginals(path + "generation_model.txt",parms);

		GenModel gen_model(parms,marginals);


		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(path + string("non_prod_alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
		sorted_alignments = read_alignments_seq_csv_score_range(path + string("non_prod_alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15
		sorted_alignments = read_alignments_seq_csv_score_range(path + string("non_prod_alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);


		gen_model.infer_model(sorted_alignments_vec , 1 , path+string("run") + string(argv[2]) + string("/"));


	}

	if(false){
		//Quick dirty inference on self gen seq
		string path (argv[1]);
		if (path[path.size()-1] != '/'){
			path+="/";
		}



		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("genomicVs_with_primers.fasta"));

		vector<pair<string,string>> d_genomic = read_genomic_fasta(path + string("genomicDs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("genomicJs_all_curated.fasta"));

		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
		d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string(argv[3]));
/*
		v_aligner.align_seqs(path + string("alignments_V.csv"), indexed_seqlist,750,true);


		d_aligner.align_seqs(path + string("alignments_D.csv") , indexed_seqlist,0,false);


		j_aligner.align_seqs(path + string("alignments_J.csv") , indexed_seqlist,80,true);*/

		//j_aligner.write_indexed_seq_csv(path + string("indexed_seq.csv") , indexed_seqlist);

		//Construct tcr model
		Gene_choice v_choice(V_gene,v_genomic);
		v_choice.set_nickname("v_choice");
		v_choice.set_priority(8);
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
		parms.add_edge(&d_choice,&j_choice);

		Model_marginals model_marginals(parms);
		model_marginals.uniform_initialize(parms);

		Single_error_rate error_rate(stod(string(argv[4])));

		parms.set_error_ratep(&error_rate);

		GenModel gen_model(parms,model_marginals);

		//Read alignments

		//vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("indexed_seq_test.csv"));
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
		//sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15
		sorted_alignments = read_alignments_seq_csv(path + string("alignments_D.csv"), D_gene , 1000 , false , indexed_seqlist , sorted_alignments);
		sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

		cout<<"num_procs: "<<omp_get_num_procs()<<endl;
		cout<<"num_threads: "<<omp_get_num_threads()<<endl;


		gen_model.infer_model(sorted_alignments_vec , 20 , path+string("run") + string(argv[2]) + string("/"),1e-40);

	}

	if(false){
		string path (argv[1]);
		if (path[path.size()-1] != '/'){
			path+="/";
		}

		//Infer misha new beta ?

		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("genomicVs_with_primers.fasta"));

		vector<pair<string,string>> d_genomic = read_genomic_fasta(path + string("genomicDs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("genomicJs_all_curated.fasta"));

		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
		d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		//vector<pair<const int, const string>> indexed_seqlist = read_txt(path + string(argv[3]));
		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string(argv[3]));

		v_aligner.align_seqs(path + string("test/alignments_V.csv"),indexed_seqlist,50,true,INT16_MIN,-115);


		d_aligner.align_seqs(path + string("test/alignments_D.csv"),indexed_seqlist,0,false);

		j_aligner.align_seqs(path + string("test/alignments_J.csv"),indexed_seqlist,10,true,45,57);

		//j_aligner.write_indexed_seq_csv(path + string("test/indexed_seq.csv") , indexed_seqlist);

	/*	//Construct tcr model
		Gene_choice v_choice(V_gene,v_genomic);
		v_choice.set_nickname("v_choice");
		v_choice.set_priority(7);
		Gene_choice d_choice(D_gene,d_genomic);
		d_choice.set_nickname("d_gene");
		d_choice.set_priority(8);
		Gene_choice j_choice(J_gene,j_genomic);
		j_choice.set_nickname("j_choice");
		j_choice.set_priority(6);

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
		parms.add_edge(&d_choice,&j_choice);

		Model_marginals model_marginals(parms);
		model_marginals.uniform_initialize(parms);

		Single_error_rate error_rate(stod(string(argv[4])));

		parms.set_error_ratep(&error_rate);

		GenModel gen_model(parms,model_marginals);

		//Read alignments
		//vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("indexed_seq.csv"));
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
		sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15
		sorted_alignments = read_alignments_seq_csv_score_range(path + string("alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

		cout<<"num_procs: "<<omp_get_num_procs()<<endl;
		cout<<"num_threads: "<<omp_get_num_threads()<<endl;


		gen_model.infer_model(sorted_alignments_vec , 10 , path+string("run") + string(argv[2]) + string("/"));
*/
	}

	if(false){

		//Infer misha old twin alpha chain model
		string path (argv[1]);
		if (path[path.size()-1] != '/'){
			path+="/";
		}



		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("TRAV.fasta"));


		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("TRAJ.fasta"));


		//Construct tcr alpha model
		Gene_choice v_choice(V_gene,v_genomic);
		v_choice.set_nickname("v_choice");
		v_choice.set_priority(7);

		Gene_choice j_choice(J_gene,j_genomic);
		j_choice.set_nickname("j_choice");
		j_choice.set_priority(6);

		Deletion v_3_del(V_gene,Three_prime,make_pair(-4,16));//16
		v_3_del.set_nickname("v_3_del");
		v_3_del.set_priority(5);

		Deletion j_5_del(J_gene,Five_prime,make_pair(-4,18));
		j_5_del.set_nickname("j_5_del");
		j_5_del.set_priority(5);

		Insertion vj_ins(VJ_genes,make_pair(0,40));
		vj_ins.set_nickname("vj_ins");
		vj_ins.set_priority(4);


		Dinucl_markov markov_model_vj(VJ_genes);
		markov_model_vj.set_nickname("vj_dinucl");
		markov_model_vj.set_priority(3);





		Model_Parms parms;


		/*parms.read_model_parms(path  + "/run_test/initial_model.txt");

		Model_marginals model_marginals(parms);
		model_marginals.txt2marginals(path  + "/run_test/initial_marginals.txt",parms);
*/



		//Add nodes to the graph
		parms.add_event(&v_choice);

		parms.add_event(&j_choice);

		parms.add_event(&v_3_del);

		parms.add_event(&j_5_del);

		parms.add_event(&vj_ins);


		parms.add_event(&markov_model_vj);



		//Add correlations
		parms.add_edge(&v_choice,&j_choice);
		parms.add_edge(&v_choice,&v_3_del);
		parms.add_edge(&j_choice,&j_5_del);


		Model_marginals model_marginals(parms);
		model_marginals.uniform_initialize(parms);

		Single_error_rate error_rate(stod(string(argv[4])));

		parms.set_error_ratep(&error_rate);

		GenModel gen_model(parms,model_marginals);

	/*	//Read alignments
		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("TwA1_A_indx.csv.csv"));

		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		unordered_map<int,forward_list<Alignment_data>> v_alignments = v_aligner.align_seqs(indexed_seqlist,20 , true, INT16_MIN , -55 );
		v_aligner.write_alignments_seq_csv(path + string("TwA1_A_V_alig_test.csv") , v_alignments);

		unordered_map<int,forward_list<Alignment_data>> j_alignments = j_aligner.align_seqs(indexed_seqlist,10, true , 58 , 75 );
		j_aligner.write_alignments_seq_csv(path + string("TwA1_A_J_alig_test.csv") , j_alignments);
*/



		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("TwA1_A_indx.csv"));
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv(path + string("TwA1_A_V_alig_test.csv"), V_gene , 20 , false , indexed_seqlist  );//40//35

		sorted_alignments = read_alignments_seq_csv(path + string("TwA1_A_J_alig_test.csv"), J_gene , 100 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);


		gen_model.infer_model(sorted_alignments_vec , 1 , path+string("run") + string(argv[2]) + string("/"));

	}

	if(false){
		string workfolder("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/MEMORY_04-AJ-M_001-025/");
		auto indexed_seqlist_all = read_indexed_csv(workfolder + "indexed_seq.csv");
		cout<<"test"<<endl;
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>>
		best_alig;
		        best_alig = read_alignments_seq_csv_score_range(workfolder + "alignments_V.csv",
		V_gene, 0, true, indexed_seqlist_all, best_alig);
		        cout<<"test"<<endl;
		        best_alig = read_alignments_seq_csv_score_range(workfolder + "alignments_J.csv",
		J_gene, 0, true, indexed_seqlist_all, best_alig);
		        ofstream outfile_ba(workfolder + "best_alignments.txt");
		        outfile_ba << "index;sequence;V gene;V score;J gene;J score;" <<
		endl;
		        for (auto& it : best_alig){
		            outfile_ba << it.first << ";" << it.second.first << ";";
		            outfile_ba << it.second.second[V_gene][0].gene_name <<  ";" <<
		it.second.second[V_gene][0].score << ";";
		            outfile_ba << it.second.second[J_gene][0].gene_name <<  ";" <<
		it.second.second[J_gene][0].score << ";";
		            outfile_ba << endl;
		        }
	}

	if(false){

        string path (argv[1]);
        if (path[path.size()-1] != '/'){
                path+="/";
        }

        //Infer BCR model
        //Align Harlan BCR sequences
/*        vector<pair<string,string>> v_genomic = read_genomic_fasta(string("/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicVs.fasta"));

        vector<pair<string,string>> d_genomic = read_genomic_fasta(string("/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicDs.fasta"));

        vector<pair<string,string>> j_genomic = read_genomic_fasta(string("/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicJs.fasta"));*/

        vector<pair<string,string>> v_genomic = read_genomic_fasta(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicVs.fasta"));

        vector<pair<string,string>> d_genomic = read_genomic_fasta(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicDs.fasta"));

        vector<pair<string,string>> j_genomic = read_genomic_fasta(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicJs.fasta"));



        //string path = "/users/marcou/PhD/bcr_harlan/memory_non_coding/alignments/MEMORY_04-AJ-M_001-025/";


        //Construct bcr model
        Gene_choice v_choice(V_gene,v_genomic);
        v_choice.set_nickname("v_choice");
        v_choice.set_priority(8);
        Gene_choice d_choice(D_gene,d_genomic);
        d_choice.set_nickname("d_gene");
        d_choice.set_priority(5);
        Gene_choice j_choice(J_gene,j_genomic);
        j_choice.set_nickname("j_choice");
        j_choice.set_priority(7);

        Deletion v_3_del(V_gene,Three_prime,make_pair(-6,20));//16
        v_3_del.set_nickname("v_3_del");
        v_3_del.set_priority(6);
        Deletion d_5_del(D_gene,Five_prime,make_pair(-6,37));
        d_5_del.set_nickname("d_5_del");
        d_5_del.set_priority(4);
        Deletion d_3_del(D_gene,Three_prime,make_pair(-6,37));
        d_3_del.set_nickname("d_3_del");
        d_3_del.set_priority(4);
        Deletion j_5_del(J_gene,Five_prime,make_pair(-6,20));
        j_5_del.set_nickname("j_5_del");
        j_5_del.set_priority(6);

        Insertion vd_ins(VD_genes,make_pair(0,60));
        vd_ins.set_nickname("vd_ins");
        vd_ins.set_priority(3);
        Insertion dj_ins(DJ_genes,make_pair(0,60));
        dj_ins.set_nickname("dj_ins");
        dj_ins.set_priority(3);

        Dinucl_markov markov_model_vd(VD_genes);
        markov_model_vd.set_nickname("vd_dinucl");
        markov_model_vd.set_priority(2);

        Dinucl_markov markov_model_dj(DJ_genes);
        markov_model_dj.set_nickname("dj_dinucl");
        markov_model_dj.set_priority(2);




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
        parms.add_edge(&v_choice,&d_choice);
        parms.add_edge(&v_choice,&j_choice);
        parms.add_edge(&v_choice,&v_3_del);
        parms.add_edge(&j_choice,&j_5_del);
        parms.add_edge(&d_choice,&d_3_del);
        parms.add_edge(&d_choice,&d_5_del);
        parms.add_edge(&d_5_del,&d_3_del);
        parms.add_edge(&j_choice,&d_choice);

        Model_marginals model_marginals(parms);
        model_marginals.uniform_initialize(parms);

        Single_error_rate error_rate(stod(string(argv[4])));

        parms.set_error_ratep(&error_rate);

        GenModel gen_model(parms,model_marginals);

        //Read alignments
        vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("AJ_Naive_noncoding_indexed_seq.csv"));
        vector<pair<const int, const string>> indexed_seqlist_sample = sample_indexed_seq(indexed_seqlist,20000);
        cout<<"indexed_seq_list sample size:"<<indexed_seqlist_sample.size()<<endl;

        unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(path + string("AJ_Naive_noncoding_alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
        sorted_alignments = read_alignments_seq_csv(path + string("AJ_Naive_noncoding_alignments_D.csv"), D_gene , 50 , false , indexed_seqlist , sorted_alignments);//30//15
        sorted_alignments = read_alignments_seq_csv_score_range(path + string("AJ_Naive_noncoding_alignments_J.csv"), J_gene , 20 , false , indexed_seqlist , sorted_alignments);//30//20

        vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

        cout<<"num_procs: "<<omp_get_num_procs()<<endl;
        cout<<"num_threads: "<<omp_get_num_threads()<<endl;


        gen_model.infer_model(sorted_alignments_vec , 20 , path+string("run") + string(argv[2]) + string("/"),1e-35,0.001);

	}

	if(true){

		//Generate sequences from the BCR naive model

		string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/NAIVE_03-AJ-N_A_026-050/");

		Model_Parms parms;
		parms.read_model_parms(path  + "/run_dchoice_bias/iteration_15_parms.txt");

		Model_marginals marginals(parms);
		marginals.txt2marginals(path  + "/run_dchoice_bias/iteration_15.txt",parms);


		GenModel gen_model(parms,marginals);


		gen_model.generate_sequences(4000000,false,path  + "/generated_sequences/" + "AJ_naive_4000000_genseqs_noerr.csv",path  + "/generated_sequences/"  + "AJ_naive_4000000_genseqs_noerr_realizations.csv");

		/*//Cut the sequences
		for(auto seq = generated_seq.begin() ; seq != generated_seq.end() ; seq++){
			string tmp = (*seq).first.substr((*seq).first.size()-101,100);
			(*seq).first = tmp;
		}*/

		marginals.write2txt(path  + "/generated_sequences/" + "AJ_naive_generation_marginals.txt",parms);
		parms.write_model_parms(path + "/generated_sequences/"  + "AJ_naive_generation_model.txt");

		//gen_model.write_seq_real2txt(path + batch + "/gen_seq_werr/" + batch + "_1000000_genseqs_werr.csv",path + batch + "/gen_seq_werr/" + batch + "_1000000_genseqs_werr_realizations.csv",generated_seq);

	}

	if(false){
		vector<pair<string,string>> v_genomic = read_genomic_fasta(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicVs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/GEN_DATA/AJ_alleles_final/genomicJs.fasta"));

		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-4,-4,-4,1,-4,-4,1,-4,1,-4,-1,-1,-1,-2,-4,5,-4,-4,-4,1,-4,1,1,-4,-1,-4,-1,-1,-2,-4,-4,5,-4,1,-4,1,-4,1,-4,-1,-1,-4,-1,-2,-4,-4,-4,5,-4,1,1,-4,-4,1,-1,-1,-1,-4,-2,1,-4,1,-4,-1,-4,-2,-2,-2,-2,-3,-1,-3,-1,-1,-4,1,-4,1,-4,-1,-2,-2,-2,-2,-1,-3,-1,-3,-1,-4,-4,1,1,-2,-2,-1,-4,-2,-2,-1,-1,-3,-3,-1,1,1,-4,-4,-2,-2,-4,-1,-2,-2,-3,-3,-1,-1,-1,-4,1,1,-4,-2,-2,-2,-2,-1,-4,-1,-3,-3,-1,-1,1,-4,-4,1,-2,-2,-2,-2,-4,-1,-3,-1,-1,-3,-1,-4,-1,-1,-1,-3,-1,-1,-3,-1,-3,-1,-2,-2,-2,-1,-1,-4,-1,-1,-1,-3,-1,-3,-3,-1,-2,-1,-2,-2,-1,-1,-1,-4,-1,-3,-1,-3,-1,-3,-1,-2,-2,-1,-2,-1,-1,-1,-1,-4,-1,-3,-3,-1,-1,-3,-2,-2,-2,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		Matrix<double> nuc44_sub_matrix(15,15,nuc44_vect);

		string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/bcr_harlan/memory_non_coding/alignments/MEMORY_04-AJ-M_001-025/test_align/");

		//V align test
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);
		vector<pair<const int, const string>> indexed_v_seq = read_fasta(path + "v_del_and_ins.fasta");

		v_aligner.align_seqs(path + string("alignments_V_test.csv"),indexed_v_seq,50,true);
		write_indexed_seq_csv(path + "v_ins_and_dels_indexed.txt",indexed_v_seq);


		//J align test
		Aligner j_aligner = Aligner(nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);
		vector<pair<const int, const string>> indexed_j_seq = read_fasta(path + "j_del_and_ins.fasta");

		j_aligner.align_seqs(path + string("alignments_J_test.csv"),indexed_j_seq,50,true);
		write_indexed_seq_csv(path + "j_ins_and_dels_indexed.txt",indexed_j_seq);
	}

	if(false){
		//Crash test the update functionality of events
		//////////////////////////////////////////////////
		string path ("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/old_twins_alpha/TwA1/test_event_update/");
		//Load an old twin alpha chain model with modified insertion distribution (only one insertion)

		Model_Parms parms;
		parms.read_model_parms(path + "generation_parms.txt");

		Model_marginals marginals(parms);
		marginals.txt2marginals(path + "generation_marginals.txt",parms);

/*		//Load an old twin alpha chain model with normal insertion distribution
		Model_Parms parms;
		parms.read_model_parms(path + "initial_parms.txt");

		Model_marginals marginals(parms);
		marginals.txt2marginals(path + "initial_marginals.txt",parms);*/

		//Generate some sequences
		GenModel gen_model(parms,marginals);
		gen_model.generate_sequences(2000,false,path +"1000_sample_werr.csv",path + "1000_sample_werr_realizations.csv");

		//Align them
        vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("../TRAV.fasta"));
        vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("../TRAJ.fasta"));

        vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("1000_sample_werr.csv"));

        double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
        Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

        Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
        v_aligner.set_genomic_sequences(v_genomic);

        Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
        j_aligner.set_genomic_sequences(j_genomic);

        v_aligner.align_seqs(path + string("1000_sample_werr_V_aligns.csv"),indexed_seqlist,20 , true, INT16_MIN , -55 );

        j_aligner.align_seqs(path + string("1000_sample_werr_J_aligns.csv"),indexed_seqlist,10, true , 58 , 75 );

        //Read the alignments
        unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv(path + string("1000_sample_werr_V_aligns.csv"), V_gene , 20 , false , indexed_seqlist  );//40//35
        sorted_alignments = read_alignments_seq_csv(path + string("1000_sample_werr_J_aligns.csv"), J_gene , 100 , false , indexed_seqlist , sorted_alignments);//30//20
        vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);


		//Re-infer model with all parameters free
        gen_model.infer_model(sorted_alignments_vec , 20 , path+string("run_allfree/"),1e-35,0.001);


		//Re infer the model with only the insertion distribution free
        auto event_map = parms.get_events_map();
        for(auto iter = event_map.begin() ; iter != event_map.end() ; iter++){
        	if((*iter).first != make_tuple(Insertion_t,VJ_genes,Undefined_side)){
        		(*iter).second->fix(true);
        		if((*iter).second->is_fixed()){
        			cout<<(*iter).second->get_name()<<endl;
        		}
        	}
        }
        parms.get_err_rate_p()->update_value(false);
        GenModel gen_model2(parms,marginals);
        gen_model2.infer_model(sorted_alignments_vec , 20 , path+string("run_insfree/"),1e-35,0.001);
	}

	if(false){
		//Infer the insertion distribution for Misha beta chains ranks 1 3 11 (DL naive and memory, cordblood)

		string path (argv[1]);
			if (path[path.size()-1] != '/'){
				path+="/";
			}
		string batch(argv[2]);

		size_t rank_size = 3000;

/////////////////////////////////////////////////////////
		//Alignment of sequences if needed
		vector<pair<const int, const string>> indx_read_from_txt = read_txt(path + batch + "_read.txt");

		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("TRBV.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("TRBJ.fasta"));
		cout<<"Genomic sequences read"<<endl;
		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);
		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		cout<<"aligning sequences"<<endl;
		v_aligner.align_seqs(path + batch + string("_V_alig.csv"),indx_read_from_txt,550,true);

		j_aligner.align_seqs(path + batch + string("_J_alig.csv"),indx_read_from_txt,160,true,70,INT16_MAX);
		write_indexed_seq_csv(path + batch + string("_indx.csv") , indx_read_from_txt);

////////////////////////////////////////////////////////

		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + batch + string("_indx.csv"));
		size_t n_sequences = indexed_seqlist.size();
		size_t chunks = n_sequences/rank_size;

		system(&("mkdir " + path +batch + "_ranksize_" + to_string(rank_size))[0]);

		for(size_t rank = 1 ; rank != chunks+1 ; rank++){
			vector<pair<const int, const string>> sub_indexed_seqlist(indexed_seqlist.begin() + (rank-1)*rank_size,indexed_seqlist.begin() + (rank)*rank_size );
			system(&("mkdir "+path + batch + "_ranksize_" + to_string(rank_size) + "/" + batch + "_rank_"+to_string(rank))[0] );
			write_indexed_seq_csv(path + batch + "_ranksize_" + to_string(rank_size) + "/"+ batch + "_rank_"+to_string(rank) +"/"+batch + "_indx_rank_"+to_string(rank)+".csv",sub_indexed_seqlist);

			//Read the alignments
			unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv(path + batch + string("_V_alig.csv"), V_gene , 0 , false , sub_indexed_seqlist  );//40//35
			sorted_alignments = read_alignments_seq_csv(path + batch + string("_J_alig.csv"), J_gene , 0 , false , sub_indexed_seqlist , sorted_alignments);//30//20

			//Do not provide D alignments
			for(vector<pair<const int , const string>>::const_iterator seq_it = sub_indexed_seqlist.begin() ; seq_it != sub_indexed_seqlist.end() ; seq_it++){
					sorted_alignments[(*seq_it).first].second[D_gene] = vector<Alignment_data>();
			}

			vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

			Model_Parms parms;
			parms.read_model_parms(path + "Azh_iteration_20_parms.txt");

			Model_marginals marginals(parms);
			marginals.txt2marginals(path + "Azh_iteration_20.txt",parms);

			//Re infer the model with only the insertion distribution free
			auto event_map = parms.get_events_map();
			for(auto iter = event_map.begin() ; iter != event_map.end() ; iter++){
				if(((*iter).first != make_tuple(Insertion_t,DJ_genes,Undefined_side)) & ((*iter).first != make_tuple(Insertion_t,VD_genes,Undefined_side)) ){
					(*iter).second->fix(true);
					if((*iter).second->is_fixed()){
						cout<<(*iter).second->get_name()<<" has been fixed"<<endl;
					}
				}
			}
			parms.get_err_rate_p()->update_value(false);
			GenModel gen_model(parms,marginals);
			try{
				gen_model.infer_model(sorted_alignments_vec , 15 , path + batch + "_ranksize_" + to_string(rank_size) + "/"+ batch + "_rank_"+to_string(rank) +"/",1e-35,0.001);
			}
			catch(exception& except){
				cout<<"Exception caught for batch:"<<batch<<" rank "<<to_string(rank)<<endl;
				continue;
			}

	}

}

if(false){
		//Compute read likelihoods for misha 5000 Azh productive

		string path ("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/DL_and_cord/Azh_Qinv_test/");


		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("genomicVs_with_primers.fasta"));

		vector<pair<string,string>> d_genomic = read_genomic_fasta(path + string("genomicDs.fasta"));

		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("genomicJs_all_curated.fasta"));
		cout<<"Genomic sequences read"<<endl;
		//Declare substitution matrix used for alignments(nuc44 here)
		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);
		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
		d_aligner.set_genomic_sequences(d_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		vector<pair<const int, const string>> indexed_seqlist = read_txt(path +string( "Azh_0_F1_5000_prod.txt"));
		cout<<"aligning sequences"<<endl;
		v_aligner.align_seqs(path  + string("Azh_0_F1_5000_prod_alignments_V.csv"),indexed_seqlist,50,true,INT16_MIN,-115);


		d_aligner.align_seqs(path  + string("Azh_0_F1_5000_prod_alignments_D.csv"),indexed_seqlist,0,false);

		j_aligner.align_seqs(path  + string("Azh_0_F1_5000_prod_alignments_J.csv"),indexed_seqlist,10,true,45,57);
		write_indexed_seq_csv(path  + string("Azh_0_F1_5000_prod_indexed_seq.csv") , indexed_seqlist);

		Model_Parms parms;
		cout<<"Reading model parms..."<<endl;
		parms.read_model_parms(path  + "../Azh_iteration_20_parms.txt");
		Model_marginals marginals(parms);
		cout<<"Reading marginals..."<<endl;
		marginals.txt2marginals(path  + "../Azh_iteration_20.txt",parms);



		GenModel gen_model(parms,marginals);


		cout<<"Reading alignments..."<<endl;
		//vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string(argv[3]) + string("_nonprod_indexed_seq.csv"));

		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(path  + string("Azh_0_F1_5000_prod_alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35


		sorted_alignments = read_alignments_seq_csv_score_range(path  + string("Azh_0_F1_5000_prod_alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15

		sorted_alignments = read_alignments_seq_csv_score_range(path  + string("Azh_0_F1_5000_prod_alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20



		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);
		cout<<"Computing probabilities..."<<endl;
		gen_model.infer_model(sorted_alignments_vec , 1 , path+ string("/"),1e-35,1e-3);




}

if(false){

		//Infer misha new twin alpha chain model
		string path (argv[1]);
		if (path[path.size()-1] != '/'){
			path+="/";
		}
		string batch(argv[2]);



		vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("TRAV.fasta"));


		vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("TRAJ.fasta"));


		//Construct tcr alpha model
		Gene_choice v_choice(V_gene,v_genomic);
		v_choice.set_nickname("v_choice");
		v_choice.set_priority(7);

		Gene_choice j_choice(J_gene,j_genomic);
		j_choice.set_nickname("j_choice");
		j_choice.set_priority(6);

		Deletion v_3_del(V_gene,Three_prime,make_pair(0,16));//16
		v_3_del.set_nickname("v_3_del");
		v_3_del.set_priority(5);

		Deletion j_5_del(J_gene,Five_prime,make_pair(0,18));
		j_5_del.set_nickname("j_5_del");
		j_5_del.set_priority(5);

		Insertion vj_ins(VJ_genes,make_pair(0,40));
		vj_ins.set_nickname("vj_ins");
		vj_ins.set_priority(4);


		Dinucl_markov markov_model_vj(VJ_genes);
		markov_model_vj.set_nickname("vj_dinucl");
		markov_model_vj.set_priority(3);
		markov_model_vj.fix(true);





		Model_Parms parms;


		/*parms.read_model_parms(path  + "/run_test/initial_model.txt");

		Model_marginals model_marginals(parms);
		model_marginals.txt2marginals(path  + "/run_test/initial_marginals.txt",parms);
*/



		//Add nodes to the graph
		parms.add_event(&v_choice);

		parms.add_event(&j_choice);

		parms.add_event(&v_3_del);

		parms.add_event(&j_5_del);

		parms.add_event(&vj_ins);


		parms.add_event(&markov_model_vj);



		//Add correlations
		parms.add_edge(&v_choice,&j_choice);
		parms.add_edge(&v_choice,&v_3_del);
		parms.add_edge(&j_choice,&j_5_del);


		Model_marginals model_marginals(parms);
		model_marginals.uniform_initialize(parms);

		Single_error_rate error_rate(stod(string(argv[4])));

		parms.set_error_ratep(&error_rate);

		GenModel gen_model(parms,model_marginals);

	/*	//Read alignments
		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + string("TwA1_A_indx.csv.csv"));

		double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
		Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);

		//Instantiate aligner with substitution matrix and gap penalty
		Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
		v_aligner.set_genomic_sequences(v_genomic);

		Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
		j_aligner.set_genomic_sequences(j_genomic);

		unordered_map<int,forward_list<Alignment_data>> v_alignments = v_aligner.align_seqs(indexed_seqlist,20 , true, INT16_MIN , -55 );
		v_aligner.write_alignments_seq_csv(path + string("TwA1_A_V_alig_test.csv") , v_alignments);

		unordered_map<int,forward_list<Alignment_data>> j_alignments = j_aligner.align_seqs(indexed_seqlist,10, true , 58 , 75 );
		j_aligner.write_alignments_seq_csv(path + string("TwA1_A_J_alig_test.csv") , j_alignments);
*/



		vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path + batch + "/"+batch+"_indx.csv");
		unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv(path  + batch +  "/" + batch + "_V_alig.csv", V_gene , 0 , false , indexed_seqlist  );//40//35

		sorted_alignments = read_alignments_seq_csv(path  + batch +  "/" + batch + "_J_alig.csv" , J_gene , 0 , false , indexed_seqlist , sorted_alignments);//30//20

		vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);


		gen_model.infer_model(sorted_alignments_vec , 20 , path+ batch + string("/run_") + string(argv[3]) + string("/") , 1e-30);

	}

	if(false){

		//string batch("TV_F");
		string batch(argv[1]);
		cout<<"processing "<<batch<<endl;
/*		string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/old_twins_alpha/");

		Model_Parms parms;
		parms.read_model_parms(path + batch + "/run_lowprobthresh/iteration_20_parms.txt");

		Model_marginals marginals(parms);
		marginals.txt2marginals(path + batch + "/run_lowprobthresh/iteration_20.txt",parms);*/

		string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/new_twins_beta/");

		Model_Parms parms;
		parms.read_model_parms(path + batch + "/run_clust/iteration_20_parms.txt");

		Model_marginals marginals(parms);
		marginals.txt2marginals(path + batch + "/run_clust/iteration_20.txt",parms);


		GenModel gen_model(parms,marginals);


		gen_model.generate_sequences(4000000,false,path + batch + "/gen_seq_noerr/" + batch + "_4000000_genseqs_noerr.csv",path + batch + "/gen_seq_noerr/" + batch + "_4000000_genseqs_noerr_realizations.csv");

		/*//Cut the sequences
		for(auto seq = generated_seq.begin() ; seq != generated_seq.end() ; seq++){
			string tmp = (*seq).first.substr((*seq).first.size()-101,100);
			(*seq).first = tmp;
		}*/

		marginals.write2txt(path + batch + "/gen_seq_noerr/" + batch + "_generation_marginals.txt",parms);
		parms.write_model_parms(path + batch + "/gen_seq_noerr/" + batch + "_generation_model.txt");

		/*//Flatten dinucleotide markov model
		auto event_map = parms.get_events_map();
		for(auto iter = event_map.begin() ; iter != event_map.end() ; iter++){
			if(((*iter).first == make_tuple(Dinuclmarkov_t,VJ_genes,Undefined_side)) ){
				marginals.flatten((*iter).second,parms);
				if((*iter).second->is_fixed()){
					cout<<(*iter).second->get_name()<<" has been flatten"<<endl;
				}
			}
		}
		GenModel genmodel_bis(parms,marginals);
		genmodel_bis.generate_sequences(4000000,false,path + batch + "/gen_seq_noerr/" + batch + "_flat_4000000_genseqs_noerr.csv",path + batch + "/gen_seq_noerr/" + batch + "_flat_4000000_genseqs_noerr_realizations.csv");
		marginals.write2txt(path + batch + "/gen_seq_noerr/" + batch + "flat_generation_marginals.txt",parms);
		parms.write_model_parms(path + batch + "/gen_seq_noerr/" + batch + "flat_generation_model.txt");*/

	}

	if(false){

			//string batch("TV_F");
			string batch(argv[1]);
			cout<<"processing "<<batch<<endl;
	/*		string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/old_twins_alpha/");

			Model_Parms parms;
			parms.read_model_parms(path + batch + "/run_lowprobthresh/iteration_20_parms.txt");

			Model_marginals marginals(parms);
			marginals.txt2marginals(path + batch + "/run_lowprobthresh/iteration_20.txt",parms);*/

			string path("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/misha_twins/new_twins_alpha/");

			Model_Parms parms;
			parms.read_model_parms(path  + "/yuval_parms_test/"+batch+"_Q_parms.txt");

			Model_marginals marginals(parms);
			marginals.txt2marginals(path  + "/yuval_parms_test/"+batch+"_Q.txt",parms);


			GenModel gen_model(parms,marginals);


			gen_model.generate_sequences(4000000,false,path + batch + "/gen_seq_noerr/" + batch + "_yuval_4000000_genseqs_noerr.csv",path + batch + "/gen_seq_noerr/" + batch + "_yuval_4000000_genseqs_noerr_realizations.csv");

			/*//Cut the sequences
			for(auto seq = generated_seq.begin() ; seq != generated_seq.end() ; seq++){
				string tmp = (*seq).first.substr((*seq).first.size()-101,100);
				(*seq).first = tmp;
			}*/

			marginals.write2txt(path + batch + "/gen_seq_noerr/" + batch + "yuval_generation_marginals.txt",parms);
			parms.write_model_parms(path + batch + "/gen_seq_noerr/" + batch + "yuval_generation_model.txt");

			/*//Flatten dinucleotide markov model
			auto event_map = parms.get_events_map();
			for(auto iter = event_map.begin() ; iter != event_map.end() ; iter++){
				if(((*iter).first == make_tuple(Dinuclmarkov_t,VJ_genes,Undefined_side)) ){
					marginals.flatten((*iter).second,parms);
					if((*iter).second->is_fixed()){
						cout<<(*iter).second->get_name()<<" has been flatten"<<endl;
					}
				}
			}
			GenModel genmodel_bis(parms,marginals);
			genmodel_bis.generate_sequences(4000000,false,path + batch + "/gen_seq_noerr/" + batch + "_flat_4000000_genseqs_noerr.csv",path + batch + "/gen_seq_noerr/" + batch + "_flat_4000000_genseqs_noerr_realizations.csv");
			marginals.write2txt(path + batch + "/gen_seq_noerr/" + batch + "flat_generation_marginals.txt",parms);
			parms.write_model_parms(path + batch + "/gen_seq_noerr/" + batch + "flat_generation_model.txt");*/

		}

	if(false){
		//Infer the insertion distribution for Misha beta chains on Aging data

			string path (argv[1]);
				if (path[path.size()-1] != '/'){
					path+="/";
				}
			string batch(argv[2]);


			/////////////////////////////////////////////////////////
			//Alignment of sequences if needed
			vector<pair<const int, const string>> indx_read_from_txt = read_txt(path +"processed/" + batch + "/" + batch + "_reads.txt");

			vector<pair<string,string>> v_genomic = read_genomic_fasta(path + string("TRBV.fasta"));

			vector<pair<string,string>> j_genomic = read_genomic_fasta(path + string("TRBJ.fasta"));
			cout<<"Genomic sequences read"<<endl;
			//Declare substitution matrix used for alignments(nuc44 here)
			double nuc44_vect [] = {5,-14,-14,-14 , -14 ,5,-14,-14 , -14,-14,5,-14 , -14,-14,-14,5};
			Matrix<double> nuc44_sub_matrix(4,4,nuc44_vect);
			//Instantiate aligner with substitution matrix and gap penalty
			Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
			v_aligner.set_genomic_sequences(v_genomic);

			Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
			j_aligner.set_genomic_sequences(j_genomic);

			cout<<"aligning sequences"<<endl;
			v_aligner.align_seqs(path +"processed/" + batch + "/" + batch + string("_V_alig.csv"),indx_read_from_txt,550,true);

			j_aligner.align_seqs(path +"processed/" + batch + "/" + batch + string("_J_alig.csv"),indx_read_from_txt,160,true,70,INT16_MAX);
			write_indexed_seq_csv(path +"processed/" + batch + "/" + batch + string("_indx.csv") , indx_read_from_txt);

	////////////////////////////////////////////////////////



			//Read the alignments
			unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv(path +"processed/" + batch + "/" + batch + string("_V_alig.csv"), V_gene , 0 , false , indx_read_from_txt  );//40//35
			sorted_alignments = read_alignments_seq_csv(path +"processed/" + batch + "/" + batch + string("_J_alig.csv"), J_gene , 0 , false , indx_read_from_txt , sorted_alignments);//30//20

			//Do not provide D alignments
			for(vector<pair<const int , const string>>::const_iterator seq_it = indx_read_from_txt.begin() ; seq_it != indx_read_from_txt.end() ; seq_it++){
					sorted_alignments[(*seq_it).first].second[D_gene] = vector<Alignment_data>();
			}

			vector<pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

			Model_Parms parms;
			parms.read_model_parms(path + "Azh_iteration_20_parms.txt");

			Model_marginals marginals(parms);
			marginals.txt2marginals(path + "Azh_iteration_20.txt",parms);

			//Re infer the model with only the insertion distribution free
			auto event_map = parms.get_events_map();
			for(auto iter = event_map.begin() ; iter != event_map.end() ; iter++){
				if(((*iter).first != make_tuple(Insertion_t,DJ_genes,Undefined_side)) & ((*iter).first != make_tuple(Insertion_t,VD_genes,Undefined_side)) ){
					(*iter).second->fix(true);
					if((*iter).second->is_fixed()){
						cout<<(*iter).second->get_name()<<" has been fixed"<<endl;
					}
				}
			}
			parms.get_err_rate_p()->update_value(false);
			GenModel gen_model(parms,marginals);
			try{
				gen_model.infer_model(sorted_alignments_vec , 15 , path + "processed/" + batch +"/" + string(argv[3]) +"/" ,1e-35,0.001);
			}
			catch(exception& except){
				cout<<"Exception caught for batch:"<<batch<<endl;
			}


	}





	return 0;


}



