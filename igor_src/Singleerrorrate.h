/*
 * Singleerrorrate.h
 *
 *  Created on: Jan 23, 2015
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

#ifndef SINGLEERRORRATE_H_
#define SINGLEERRORRATE_H_

#include "Errorrate.h"
#include "Utils.h"
#include <math.h>

//Debug
#include <iostream>

/**
 * \class Single_error_rate Singleerrorrate.h
 * \brief Independent single nucleotide error model.
 * \author Q.Marcou
 * \version 1.0
 *
 * Simplest instance of the ErrorRate family. Models errors/mutations as a Bernouilli process with a global rate independent of position and context.
 */
class Single_error_rate: public Error_rate {
public:
	Single_error_rate();
	Single_error_rate(double);
	virtual ~Single_error_rate();
	double compare_sequences_error_prob( double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& , Mismatch_vectors_map& , double& , double& );
	void update();
	void add_to_norm_counter();
	void clean_seq_counters();
	Single_error_rate operator+(Single_error_rate);
	Single_error_rate& operator+=(Single_error_rate);
	void write2txt(std::ofstream&);
	std::shared_ptr<Error_rate> copy()const;
	std::string type() const {return "SingleErrorRate";}
	Error_rate* add_checked (Error_rate*);
	const double& get_err_rate_upper_bound(size_t,size_t) ;
	void build_upper_bound_matrix(size_t,size_t);
	int get_number_non_zero_likelihood_seqs() const{return number_seq;};
	std::queue<int>  generate_errors(std::string& , std::default_random_engine&) const;






private:

	double model_rate;
	double normalized_counter;
	long double seq_weighted_er;
	int number_errors;
	int genomic_nucl;
	long double temp2;
	long double temp;


	int subseq_compare_err_num(const std::string& , const std::string&);
	//TODO use seq likelihood to extract the likelihood of the model on the fly


};

#endif /* SINGLEERRORRATE_H_ */
