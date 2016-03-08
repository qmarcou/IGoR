/*
 * Hypermutationglobalerrorrate.cpp
 *
 *  Created on: Mar 8, 2016
 *      Author: quentin
 */

#include "Hypermutationglobalerrorrate.h"

using namespace std;

Hypermutation_global_errorrate::Hypermutation_global_errorrate() {
	// TODO Auto-generated constructor stub

}

Hypermutation_global_errorrate::~Hypermutation_global_errorrate() {
	// TODO Auto-generated destructor stub
}

Error_rate* Hypermutation_global_errorrate::copy()const{

}

Error_rate* Hypermutation_global_errorrate::add_checked(Error_rate* err_r){

}

double Hypermutation_global_errorrate::get_err_rate_upper_bound() const{

}

double Hypermutation_global_errorrate::compare_sequences_error_prob (double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& events_map , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){

}

queue<int> Hypermutation_global_errorrate::generate_errors(string& generated_seq , default_random_engine& generator) const{

}

void Hypermutation_global_errorrate::update(){

}

void Hypermutation_global_errorrate::add_to_norm_counter(){

}

void Hypermutation_global_errorrate::clean_seq_counters(){

}

void Hypermutation_global_errorrate::write2txt(ofstream& outfile){

}
