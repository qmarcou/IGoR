/*
 * Pgencounter.h
 *
 *  Created on: Aug 19, 2016
 *      Author: quentin
 */

#ifndef GENERATIVE_MODEL_SRC_PGENCOUNTER_H_
#define GENERATIVE_MODEL_SRC_PGENCOUNTER_H_

#include "Counter.h"
#include <unordered_map>

class Pgen_counter: public Counter {
public:
	Pgen_counter();
	virtual ~Pgen_counter();

	void initialize_counter(const Model_Parms& , const Model_marginals&) ;

	void count_scenario(double , double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& );

	void dump_sequence_data(int , int);


private:
	std::ofstream output_pgen_file;
	std::unordered_map<std::string,std::pair<double,double>> sequence_Pgens_map;
	std::string scenario_resulting_sequence;

	bool output_sequences;

	double read_likelihood;

	bool v_gene;
	bool d_gene;
	bool j_gene;
	bool vd_ins;
	bool dj_ins;
	bool vj_ins;
};

#endif /* GENERATIVE_MODEL_SRC_PGENCOUNTER_H_ */
