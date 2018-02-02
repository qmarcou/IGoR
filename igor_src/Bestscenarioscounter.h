/*
 * Bestscenarioscounter.h
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
 */

#ifndef GENERATIVE_MODEL_SRC_BESTSCENARIOSCOUNTER_H_
#define GENERATIVE_MODEL_SRC_BESTSCENARIOSCOUNTER_H_

#include "Counter.h"
#include <queue>
#include <vector>
#include <list>

#include "Genechoice.h"
#include "Deletion.h"
#include "Insertion.h"


/**
 * \class Best_scenarios_counter Bestscenarioscounter.h
 * \brief Records the N best scenarios realizations and mismatches.
 * \author Q.Marcou
 * \version 1.0
 *
 * Implementation of the Counter abstract class.
 * Records the N most likely scenario realizations and mismatches and append it to a semicolon separated file.
 */
class Best_scenarios_counter: public Counter {
public:
	Best_scenarios_counter();
	Best_scenarios_counter(size_t);
	Best_scenarios_counter(size_t , bool);
	Best_scenarios_counter(size_t , std::string);
	Best_scenarios_counter(size_t , std::string , bool);
	virtual ~Best_scenarios_counter();

	std::string type() const{return "BestScenarioCounter";};//TODO return an enum

	void initialize_counter(const Model_Parms& , const Model_marginals&) ;

	void count_scenario(long double , double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& );

	void count_sequence(double , const Model_marginals& ,const Model_Parms&);

	void add_checked(std::shared_ptr<Counter>);

	void dump_sequence_data(int , int);

	std::shared_ptr<Counter> copy() const;

	size_t n_scenarios_counted;
	std::shared_ptr<std::ofstream> output_scenario_file_ptr;

	std::queue< std::vector<int>> single_scenario_realizations_queue;

	std::list<int> single_scenario_mismatches_list;

	std::vector<std::tuple<double,std::queue<std::vector<int>>,std::list<int>>> best_scenarios_vec;

	std::forward_list<std::shared_ptr<const Rec_Event>> event_fw_list;

};


#endif /* GENERATIVE_MODEL_SRC_BESTSCENARIOSCOUNTER_H_ */
