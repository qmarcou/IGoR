/*
 * Hypermutationglobalerrorrate.h
 *
 *  Created on: Mar 8, 2016
 *      Author: quentin
 */

#ifndef HYPERMUTATIONGLOBALERRORRATE_H_
#define HYPERMUTATIONGLOBALERRORRATE_H_

#include "Errorrate.h"
#include <algorithm>
#include <array>

class Hypermutation_global_errorrate: public Error_rate {
	friend Gene_choice; //Grant friendship to access current gene realization and offset
	friend Deletion;//Grant friendship to access the current number of deletion

public:
	Hypermutation_global_errorrate(size_t,Gene_class,Gene_class,double);
	//Hypermutation_global_errorrate(size_t,Gene_class,Gene_class, ??); Constructor to read or copy the error rate
	virtual ~Hypermutation_global_errorrate();
	double compare_sequences_error_prob( double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>&  , Mismatch_vectors_map& , double& , double& );
	void update();
	void initialize(const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>&);
	void add_to_norm_counter();
	void clean_seq_counters();
	void write2txt(std::ofstream&);
	Error_rate* copy()const;
	std::string type() const {return "HypermutationGlobalErrorRate";}
	Error_rate* add_checked (Error_rate*);
	double get_err_rate_upper_bound() const;
	int get_number_non_zero_likelihood_seqs() const{return number_seq;};
	std::queue<int>  generate_errors(std::string& , std::default_random_engine&) const;

private:
	int number_seq; //FIXME check if need to remove this from here and from singleerrorrate and transfer it to errorrate

	Gene_class learn_on;
	Gene_class apply_to;
	size_t mutation_Nmer_size;
	double* ei_nucleotide_contributions;
	double Z;
	std::map<int,double> Nmer_background_proba;

	std::map<int,double> Nmer_background_proba_count;
	std::map<int,double> Nmer_SHM_proba_count;

	//Normalized coverage and error counters
	//Use C arrays (faster than maps)
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

	int*& vgene_offset_p;
	int*& dgene_offset_p;
	int*& jgene_offset_p;

	int*& vgene_real_index_p;
	int*& dgene_real_index_p;
	int*& jgene_real_index_p;

	//Get deletion values
	//TODO need to change this in order to handle multiple models
	int* v_3_del_value_p;
	int* d_5_del_value_p;
	int* d_3_del_value_p;
	int* j_5_del_value_p;
	int no_del_buffer = 0; //buffer used in case of no deletion event

	//Utility variables
	int i;//iteration utility
	int v_3_del_value_corr;//Corrected value for deletion numbers to avoid taking into account negative deletions
	int d_5_del_value_corr;
	int d_3_del_value_corr;
	int j_5_del_value_corr;



};

#endif /* HYPERMUTATIONGLOBALERRORRATE_H_ */
