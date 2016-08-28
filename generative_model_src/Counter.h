/*
 * Counter.h
 *
 *  Created on: Aug 19, 2016
 *      Author: quentin
 */

#ifndef GENERATIVE_MODEL_SRC_COUNTER_H_
#define GENERATIVE_MODEL_SRC_COUNTER_H_

#include <string>
#include <memory>
#include <fstream>
#include "Model_marginals.h"
#include "Model_Parms.h"
#include "Rec_Event.h"
#include "IntStr.h"



class Counter {
public:
	Counter();
	Counter(std::string);
	virtual ~Counter();

	virtual std::string type() const =0;

	virtual void initialize_counter(const Model_Parms& , const Model_marginals&) = 0;

	virtual void count_scenario(double , double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& );

	virtual void count_sequence();


	virtual void add_to_counter(std::shared_ptr<Counter>);
	virtual std::shared_ptr<Counter> add_checked(std::shared_ptr<Counter>) =0;

	virtual void dump_sequence_data(int , int);
	virtual void dump_data_summary();

	bool is_last_iter_only() const {return last_iter_only;}

	virtual std::shared_ptr<Counter> copy() const = 0;

protected:
	std::string path_to_file;
	bool last_iter_only;
	bool fstreams_created;
	//TODO create a unique identifier of the counter? Make something up to prevent to have twice the same counter??
};


#endif /* GENERATIVE_MODEL_SRC_COUNTER_H_ */
