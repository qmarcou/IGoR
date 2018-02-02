/*
 * Errorrate.h
 *
 *  Created on: Jan 22, 2015
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

#ifndef ERRORRATE_H_
#define ERRORRATE_H_

#include "Utils.h"
#include "IntStr.h"
#include <unordered_map>
#include <utility>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <random>
#include <queue>
#include <memory>

//Debug
#include <iostream>
#include <cmath>

//Forward declare Rec_event
class Rec_Event;

/**
 * \class Error_rate Error_rate.h
 * \brief Abstract class for generic error models behavior.
 * \author Q.Marcou
 * \version 1.0
 *
 * Base class for defining different error models such as additive or non-additive hypermutation models.
 * Errors are assessed when all RecEvent iterate have been processed (terminal leaf of the scenario tree)
 *
 */
class Error_rate {
public:
	Error_rate();
	virtual ~Error_rate();
	virtual double compare_sequences_error_prob( double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& , double& , double& )=0;
	virtual void update()=0;
	virtual void initialize(const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&);
	bool is_updated() const {return updated;}
	void update_value(bool update_status) {updated = update_status;};
	virtual void add_to_norm_counter()=0;
	virtual void clean_seq_counters()=0;
	void norm_weights_by_seq_likelihood(Marginal_array_p&, const size_t, const double seq_weight=1);
	virtual void write2txt(std::ofstream&)=0;
	virtual std::shared_ptr<Error_rate> copy() const = 0;
	virtual std::string type() const =0;
	virtual Error_rate* add_checked(Error_rate*) = 0;
	double get_model_likelihood() const{return model_log_likelihood;}
	double get_seq_likelihood() const{return seq_likelihood;}
	double get_seq_probability() const{return seq_probability;}
	double get_seq_mean_error_number() const;
	virtual const double& get_err_rate_upper_bound(size_t,size_t) =0;
	virtual void build_upper_bound_matrix(size_t,size_t) =0;
	virtual int get_number_non_zero_likelihood_seqs() const =0;
	virtual std::queue<int>  generate_errors(std::string& , std::default_random_engine&) const =0;
	void set_viterbi_run(bool viterbi_like){viterbi_run = viterbi_like;}
	int debug_number_scenarios;


protected:
	bool updated;
	long double model_log_likelihood;
	int number_seq;
	long double seq_likelihood;
	double seq_mean_error_number;
	long double scenario_new_proba;//TODO rename this guy
	long double seq_probability; //Probability of generating one sequence without taking errors into account
	bool viterbi_run;
	Matrix<double> upper_bound_proba_mat; //Store the value of the error cost of i errors and j no errors
	size_t max_err;
	size_t max_noerr;

};

void add_to_err_rate(Error_rate*,Error_rate*);


#endif /* ERRORRATE_H_ */
