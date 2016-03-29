/*
 * Singleerrorrate.h
 *
 *  Created on: Jan 23, 2015
 *      Author: quentin
 */

#ifndef SINGLEERRORRATE_H_
#define SINGLEERRORRATE_H_

#include "Errorrate.h"
#include "Utils.h"
#include <math.h>

//Debug
#include <iostream>

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
	Error_rate* copy()const;
	std::string type() const {return "SingleErrorRate";}
	Error_rate* add_checked (Error_rate*);
	double get_err_rate_upper_bound() const;
	int get_number_non_zero_likelihood_seqs() const{return number_seq;};
	std::queue<int>  generate_errors(std::string& , std::default_random_engine&) const;






private:

	double model_rate;
	double normalized_counter;
	int number_seq;
	long double seq_weighted_er;
	int number_errors;
	int genomic_nucl;
	long double scenario_new_proba;
	long double temp2;
	long double temp;


	int subseq_compare_err_num(const std::string& , const std::string&);
	//TODO use seq likelihood to extract the likelihood of the model on the fly


};

#endif /* SINGLEERRORRATE_H_ */
