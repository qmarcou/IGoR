/*
 * Hypermutationglobalerrorrate.h
 *
 *  Created on: Mar 8, 2016
 *      Author: quentin
 */

#ifndef HYPERMUTATIONGLOBALERRORRATE_H_
#define HYPERMUTATIONGLOBALERRORRATE_H_

#include "Errorrate.h"
#include "Genechoice.h"
#include "Deletion.h"
#include <algorithm>
#include <array>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <memory>


class Hypermutation_global_errorrate: public Error_rate {

public:
	Hypermutation_global_errorrate(size_t,Gene_class,Gene_class,double);
	Hypermutation_global_errorrate(size_t,Gene_class,Gene_class,double,std::vector<double>);
	//Hypermutation_global_errorrate(size_t,Gene_class,Gene_class, ??); Constructor to read or copy the error rate
	virtual ~Hypermutation_global_errorrate();
	double compare_sequences_error_prob( double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& , double& , double& );
	void update();
	void initialize(const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&);
	void add_to_norm_counter();
	void clean_seq_counters();
	void clean_all_counters();
	void write2txt(std::ofstream&);
	std::shared_ptr<Error_rate> copy()const;
	std::string type() const {return "HypermutationGlobalErrorRate";}
	Hypermutation_global_errorrate& operator+=(Hypermutation_global_errorrate);
	Error_rate* add_checked (Error_rate*);
	double get_err_rate_upper_bound() const;
	int get_number_non_zero_likelihood_seqs() const{return number_seq;};
	std::queue<int>  generate_errors(std::string& , std::default_random_engine&) const;
	unsigned generate_random_contributions(double);

private:
	void update_Nmers_proba(int,int,double);
	void compute_P_SHM_and_BG();
	double compute_Nmer_unorm_score(int*);

	void introduce_uniform_transversion(char&, std::default_random_engine& , std::uniform_real_distribution<double>&) const;

	Gene_class learn_on;
	Gene_class apply_to;
	size_t mutation_Nmer_size;
	//std::unique_ptr<double[]> ei_nucleotide_contributions;
	double* ei_nucleotide_contributions;
	double R;
	//std::map<int,double> Nmer_background_proba;
	double* Nmer_mutation_proba;

	double* Nmer_P_SHM;
	double* Nmer_P_BG;

	size_t alphabet_size = 4;

	//std::map<int,double> Nmer_background_proba_count;
	//std::map<int,double> Nmer_SHM_proba_count;


	//Normalized coverage and error counters
	//# V D and J possible realizations
	std::shared_ptr<Gene_choice> v_gene_event_p;
	size_t n_v_real;
	std::unordered_map<std::string , Event_realization> v_realizations;
	std::shared_ptr<Gene_choice> d_gene_event_p;
	size_t n_d_real;
	std::unordered_map<std::string , Event_realization> d_realizations;
	std::shared_ptr<Gene_choice> j_gene_event_p;
	size_t n_j_real;
	std::unordered_map<std::string , Event_realization> j_realizations;
	//Use C arrays (faster than maps)
	//size_t is the length of the sequence
	//double* is the pointer to the array for coverage/error per nucleotide
	std::pair<size_t,double*>* v_gene_nucleotide_coverage_p;
	std::pair<size_t,double*>* v_gene_per_nucleotide_error_p;
	std::pair<size_t,double*>* d_gene_nucleotide_coverage_p;
	std::pair<size_t,double*>* d_gene_per_nucleotide_error_p;
	std::pair<size_t,double*>* j_gene_nucleotide_coverage_p;
	std::pair<size_t,double*>* j_gene_per_nucleotide_error_p;

	//One seq coverage and error counters
	std::pair<size_t,double*>* v_gene_nucleotide_coverage_seq_p;
	std::pair<size_t,double*>* v_gene_per_nucleotide_error_seq_p;
	std::pair<size_t,double*>* d_gene_nucleotide_coverage_seq_p;
	std::pair<size_t,double*>* d_gene_per_nucleotide_error_seq_p;
	std::pair<size_t,double*>* j_gene_nucleotide_coverage_seq_p;
	std::pair<size_t,double*>* j_gene_per_nucleotide_error_seq_p;


	bool apply_to_v;
	bool apply_to_d;
	bool apply_to_j;

	bool learn_on_v;
	bool learn_on_d;
	bool learn_on_j;

	bool v_gene;
	bool d_gene;
	bool j_gene;
	bool vd_ins;
	bool dj_ins;
	bool vj_ins;

	const int** vgene_offset_p;
	const int** dgene_offset_p;
	const int** jgene_offset_p;

	const int** vgene_real_index_p;
	const int** dgene_real_index_p;
	const int** jgene_real_index_p;

	//Get deletion values
	//TODO need to change this in order to handle multiple models
	const int* v_3_del_value_p;
	const int* d_5_del_value_p;
	const int* d_3_del_value_p;
	const int* j_5_del_value_p;
	const int no_del_buffer = 0; //buffer used in case of no deletion event

	//Utility speed variables
	mutable int i;//iteration utility
	mutable int j;
	int v_3_del_value_corr;//Corrected value for deletion numbers to avoid taking into account negative deletions
	int d_5_del_value_corr;
	int d_3_del_value_corr;
	int j_5_del_value_corr;
	mutable double* tmp_cov_p;
	mutable double* tmp_err_p;
	int tmp_corr_len;
	int tmp_len_util;
	double scenario_new_proba;
	std::string scenario_resulting_sequence;

	std::vector<size_t> adressing_vector;
	mutable std::queue<size_t> current_Nmer;
	size_t largest_nuc_adress;
	mutable int tmp_int_nt;
	mutable int Nmer_index;
	std::vector<int>::const_iterator current_mismatch;

	double* debug_v_seq_coverage;
	double* debug_mismatch_seq_coverage;
	std::string debug_current_string;
///////////////
	//Ddddebug shit
	double* debug_one_seq_Nmer_N_SHM;
	double* debug_one_seq_Nmer_N_bg;
	double* debug_Nmer_N_SHM;
	double* debug_Nmer_N_bg;

};

#endif /* HYPERMUTATIONGLOBALERRORRATE_H_ */
