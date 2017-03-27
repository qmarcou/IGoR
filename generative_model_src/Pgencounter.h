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
