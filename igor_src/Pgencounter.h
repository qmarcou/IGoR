/*
 * Pgencounter.h
 *
 *  Created on: Aug 19, 2016
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

#ifndef GENERATIVE_MODEL_SRC_PGENCOUNTER_H_
#define GENERATIVE_MODEL_SRC_PGENCOUNTER_H_

#include "Counter.h"
#include <unordered_map>

/**
 * \class Pgen_counter Pgencounter.h
 * \brief Estimates sequences generation probability.
 * \author Q.Marcou
 * \version 1.0
 *
 * This Counter implements an estimator for the generation probability of evaluated sequences.
 * Alternatively the counter can record the probability of generation of putative ancestor (unmutated/error free) sequences and their associated posterior probability.
 */
class Pgen_counter: public Counter {
public:
	Pgen_counter();
	Pgen_counter(std::string);
	Pgen_counter(std::string , bool , bool do_output_sequences=false);
	virtual ~Pgen_counter();

	std::string type() const{return "PgenCounter";};//TODO return an enum

	void initialize_counter(const Model_Parms& , const Model_marginals&) ;

	void count_scenario(long double , double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& );

	void dump_sequence_data(int , int);

	void add_checked(std::shared_ptr<Counter>);

	std::shared_ptr<Counter> copy() const;




private:

	bool output_sequences;
	bool output_Pgen_estimator;

	std::shared_ptr<std::ofstream> output_pgen_file_ptr;
	std::unordered_map<Int_Str,std::pair<double,long double>> sequence_Pgens_map;
	Int_Str scenario_resulting_sequence;

	long double read_likelihood;

	bool v_gene;
	bool d_gene;
	bool j_gene;
	bool vd_ins;
	bool dj_ins;
	bool vj_ins;
};

#endif /* GENERATIVE_MODEL_SRC_PGENCOUNTER_H_ */
