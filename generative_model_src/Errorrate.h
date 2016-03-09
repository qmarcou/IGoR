/*
 * Errorrate.h
 *
 *  Created on: Jan 22, 2015
 *      Author: quentin
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

class Error_rate {
public:
	Error_rate();
	virtual ~Error_rate();
	virtual double compare_sequences_error_prob( double ,const std::string& , Seq_type_str_p_map& , const Seq_offsets_map& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&  , Mismatch_vectors_map& , double& , double& )=0;
	virtual void update()=0;
	virtual void initialize(const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>&);
	bool is_updated() const {return updated;}
	void update_value(bool update_status) {updated = update_status;};
	virtual void add_to_norm_counter()=0;
	virtual void clean_seq_counters()=0;
	void norm_weights_by_seq_likelihood(Marginal_array_p, const size_t, const double seq_weight=1);
	virtual void write2txt(std::ofstream&)=0;
	virtual std::shared_ptr<Error_rate> copy() const = 0;
	virtual std::string type() const =0;
	virtual Error_rate* add_checked(Error_rate*) = 0;
	double get_model_likelihood() const{return model_log_likelihood;}
	double get_seq_likelihood() const{return seq_likelihood;}
	double get_seq_probability() const{return seq_probability;}
	double get_seq_mean_error_number() const;
	virtual double get_err_rate_upper_bound() const =0;
	virtual int get_number_non_zero_likelihood_seqs() const =0;
	virtual std::queue<int>  generate_errors(std::string& , std::default_random_engine&) const =0;
	void set_viterbi_run(bool viterbi_like){viterbi_run = viterbi_like;}
	int debug_number_scenarios;


protected:
	bool updated;
	long double model_log_likelihood;
	long double seq_likelihood;
	double seq_mean_error_number;
	long double scenario_new_proba;//TODO rename this guy
	long double seq_probability; //Probability of generating one sequence without taking errors into account
	bool viterbi_run;

};

void add_to_err_rate(Error_rate*,Error_rate*);


#endif /* ERRORRATE_H_ */
