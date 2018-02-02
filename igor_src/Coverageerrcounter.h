/*
 * Coverageerrcounter.h
 *
 *  Created on: Oct 24, 2016
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
 */

#ifndef COVERAGEERRCOUNTER_H_
#define COVERAGEERRCOUNTER_H_

#include "Counter.h"
#include <string>
/**
 * \class Coverage_err_counter Coverageerrcounter.h
 * \brief Records the number of time each genomic position is observed with or without error.
 * \author Q.Marcou
 * \version 1.0
 *
 * The Coverage_err_counter allows to record the number of times each genomic site is observed (coverage) and how many times a mismatch has been observed on it (error/mutation).
 * The recording can be made at the single position level, as well as joint over positions duet,triplet etc (e.g the number of times two nucleotides were observed in the same scenario)
 */
class Coverage_err_counter: public Counter {
public:

	Coverage_err_counter(Gene_class);
	Coverage_err_counter(Gene_class,bool,bool);
	Coverage_err_counter(std::string,Gene_class,bool);
	Coverage_err_counter(std::string,Gene_class,size_t,bool,bool);
	virtual ~Coverage_err_counter();

	std::string type() const{return "CoverageErrCounter";};//TODO return an enum

	void initialize_counter(const Model_Parms& , const Model_marginals&) ;

	void count_scenario(long double , double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& );

	void count_sequence(double , const Model_marginals& ,const Model_Parms&);

	void add_checked(std::shared_ptr<Counter>);

	void dump_sequence_data(int , int);
	void dump_data_summary(int);

	std::shared_ptr<Counter> copy() const;

private:

	void allocate_coverage_and_errors_arrays(size_t,const std::unordered_map<std::string , Event_realization>,std::pair<size_t,double*>*&,std::pair<size_t,double*>*&,std::pair<size_t,double*>*&,std::pair<size_t,double*>*&);
	void deallocate_coverage_and_errors_arrays(size_t,const std::unordered_map<std::string , Event_realization>,std::pair<size_t,double*>*&,std::pair<size_t,double*>*&,std::pair<size_t,double*>*&,std::pair<size_t,double*>*&);
	void dump_cov_and_err_arrays(int,int,std::shared_ptr<std::ofstream> , size_t , std::pair<size_t,double*>* , std::pair<size_t,double*>*);
	void normalize_and_add_cov_and_err(double& , size_t , std::pair<size_t,double*>* , std::pair<size_t,double*>* , std::pair<size_t,double*>* , std::pair<size_t,double*>* );
	void recurs_coverage_count(double scenario_seq_joint_proba , size_t N , size_t begin_bound , size_t end_bound , size_t gene_len);
	void recurs_errors_count(double scenario_seq_joint_proba , std::vector<int>& v_mismatch_list , 	const int** gene_offset_p  , size_t N , size_t begin_bound , size_t end_bound , size_t gene_len);
	void symmetrize_counter_array(double* , size_t , size_t,size_t);
	void symmetrize_counter_array_recurs(size_t , size_t , size_t* , double* , size_t);

	std::shared_ptr<std::ofstream> output_cov_err_v_file_ptr;
	std::shared_ptr<std::ofstream> output_cov_err_d_file_ptr;
	std::shared_ptr<std::ofstream> output_cov_err_j_file_ptr;

	Gene_class count_on;
	bool dump_individual_seqs;
	size_t record_Npoint_occurence;
	size_t* positions;

	//Normalized coverage and error counters
	//# V D and J possible realizations and events
		std::shared_ptr<Gene_choice> v_gene_event_p;
		size_t n_v_real;
		std::unordered_map<std::string , Event_realization> v_realizations;
		std::shared_ptr<Gene_choice> d_gene_event_p;
		size_t n_d_real;
		std::unordered_map<std::string , Event_realization> d_realizations;
		std::shared_ptr<Gene_choice> j_gene_event_p;
		size_t n_j_real;
		std::unordered_map<std::string , Event_realization> j_realizations;

	// Use C arrays (faster than maps)
	// size_t is the length of the sequence
	// double* is the pointer to the array for coverage/error per nucleotide
		std::pair<size_t,double*>* v_gene_nucleotide_coverage_p;
		std::pair<size_t,double*>* v_gene_per_nucleotide_error_p;
		std::pair<size_t,double*>* d_gene_nucleotide_coverage_p;
		std::pair<size_t,double*>* d_gene_per_nucleotide_error_p;
		std::pair<size_t,double*>* j_gene_nucleotide_coverage_p;
		std::pair<size_t,double*>* j_gene_per_nucleotide_error_p;

	// One seq coverage and error counters
		std::pair<size_t,double*>* v_gene_nucleotide_coverage_seq_p;
		std::pair<size_t,double*>* v_gene_per_nucleotide_error_seq_p;
		std::pair<size_t,double*>* d_gene_nucleotide_coverage_seq_p;
		std::pair<size_t,double*>* d_gene_per_nucleotide_error_seq_p;
		std::pair<size_t,double*>* j_gene_nucleotide_coverage_seq_p;
		std::pair<size_t,double*>* j_gene_per_nucleotide_error_seq_p;

	// Count on
		bool count_on_v;
		bool count_on_d;
		bool count_on_j;

	//Offset pointers
		const int** vgene_offset_p;
		const int** dgene_offset_p;
		const int** jgene_offset_p;

	//Realizations pointers
		const int** vgene_real_index_p;
		const int** dgene_real_index_p;
		const int** jgene_real_index_p;

	//Get deletion values
	//TODO need to change this in order to handle multiple models(?)
		const int* v_3_del_value_p;
		const int* d_5_del_value_p;
		const int* d_3_del_value_p;
		const int* j_5_del_value_p;
		const int no_del_buffer = 0; //buffer used in case of no deletion event

	//Utility
		int i;
		mutable double* tmp_cov_p;
		mutable double* tmp_err_p;
		int tmp_corr_len;
		int tmp_len_util;

};

#endif /* COVERAGEERRCOUNTER_H_ */
