/*
 * Errorscounter.h
 *
 *  Created on: Oct 5, 2017
 *      Author: quentin
 */

#ifndef IGOR_SRC_ERRORSCOUNTER_H_
#define IGOR_SRC_ERRORSCOUNTER_H_

#include "Counter.h"

class Errors_counter: public Counter {
public:
	Errors_counter();
	Errors_counter(size_t);
	Errors_counter(size_t , bool);
	Errors_counter(size_t , std::string);
	Errors_counter(size_t , std::string , bool);
	virtual ~Errors_counter();

	std::string type() const{return "ErrorsCounter";};//TODO return an enum

	void initialize_counter(const Model_Parms& , const Model_marginals&) ;

	void count_scenario(long double , double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& );

	void count_sequence(double , const Model_marginals& ,const Model_Parms&);

	void add_checked(std::shared_ptr<Counter>);

	void dump_sequence_data(int , int);

	std::shared_ptr<Counter> copy() const;

private:
	// Gene class to count errors on different genes?

	size_t n_scenarios_counted;
	bool output_scenarios;
	std::shared_ptr<std::ofstream> output_scenario_errors_file_ptr;
	std::shared_ptr<std::ofstream> output_sequence_averaged_errors_file_ptr;


	size_t scenario_n_mismatches;
	size_t scenario_n_genomic;

	double sequence_average_error_freq;
	double sequence_average_n_genomic;
	double sequence_average_n_mismatches;

	std::vector<std::tuple<double,size_t,size_t>> best_scenarios_vec;

};


#endif /* IGOR_SRC_ERRORSCOUNTER_H_ */
