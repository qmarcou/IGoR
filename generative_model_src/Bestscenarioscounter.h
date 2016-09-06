/*
 * Bestscenarioscounter.h
 *
 *  Created on: Aug 19, 2016
 *      Author: quentin
 */

#ifndef GENERATIVE_MODEL_SRC_BESTSCENARIOSCOUNTER_H_
#define GENERATIVE_MODEL_SRC_BESTSCENARIOSCOUNTER_H_

#include "Counter.h"
#include <queue>
#include <vector>
#include <forward_list>

#include "Genechoice.h"
#include "Deletion.h"
#include "Insertion.h"



class Best_scenarios_counter: public Counter {
public:
	Best_scenarios_counter();
	Best_scenarios_counter(size_t);
	virtual ~Best_scenarios_counter();

	void initialize_counter(const Model_Parms& , const Model_marginals&) ;

	void count_scenario(double , double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& );

	void count_sequence(double , const Model_marginals& ,const Model_Parms&);

	void dump_sequence_data(int , int);

	size_t n_scenarios_counted;
	std::ofstream output_scenario_file;

	std::queue< std::vector<int>> single_scenario_realizations_queue;

	std::vector<std::pair<double,std::queue<std::vector<int>>>> best_scenarios_vec;

	std::forward_list<std::shared_ptr<const Rec_Event>> event_fw_list;

};


#endif /* GENERATIVE_MODEL_SRC_BESTSCENARIOSCOUNTER_H_ */
