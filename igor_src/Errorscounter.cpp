/*
 * Errorscounter.cpp
 *
 *  Created on: Oct 5, 2017
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

#include "Errorscounter.h"


using namespace std;


Errors_counter::Errors_counter(size_t n_scenarios): Counter(), n_scenarios_counted(n_scenarios) , output_scenarios(n_scenarios_counted>0) ,
		scenario_n_mismatches(0) , scenario_n_genomic(0) , sequence_average_error_freq(0) , sequence_average_n_genomic(0) , sequence_average_n_mismatches(0){
}

Errors_counter::Errors_counter(size_t n_scenarios , bool is_last_iter_only): Errors_counter(n_scenarios) {
	this->last_iter_only = is_last_iter_only;
}

Errors_counter::Errors_counter(size_t n_scenarios , string path): Counter(path) , n_scenarios_counted(n_scenarios) , output_scenarios(n_scenarios_counted>0),
		scenario_n_mismatches(0) , scenario_n_genomic(0) , sequence_average_error_freq(0) , sequence_average_n_genomic(0) , sequence_average_n_mismatches(0){
}

/**
 * \param n_scenarios: number of scenarios to be recorded, if 0 only the average sequence mutation rate will be recorded.
 */
Errors_counter::Errors_counter( size_t n_scenarios , string path , bool is_last_iter_only) : Errors_counter(n_scenarios , path) {
	this->last_iter_only = is_last_iter_only;
}

Errors_counter::Errors_counter():Errors_counter(0,true) {
	// TODO Auto-generated constructor stub

}

Errors_counter::~Errors_counter() {
	// TODO Auto-generated destructor stub
}

void Errors_counter::count_scenario(long double scenario_seq_joint_proba , double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists){

	// Get the number of genomic nucleotides
	if(constructed_sequences.exist(V_gene_seq)){
		const vector<int>& v_seq = *constructed_sequences.at(V_gene_seq);
		this->scenario_n_genomic+=v_seq.size();
	}
	if(constructed_sequences.exist(D_gene_seq)){
		const vector<int>& d_seq = *constructed_sequences[D_gene_seq];
		this->scenario_n_genomic+=d_seq.size();
	}
	if(constructed_sequences.exist(J_gene_seq)){
		const vector<int>& j_seq = *constructed_sequences.at(J_gene_seq);
		this->scenario_n_genomic+=j_seq.size();
	}

	// Get the number of mismatches and add them to the mismatch counter
	if(mismatches_lists.exist(V_gene_seq)){
		const vector<int>& v_mismatch_list = *mismatches_lists.at(V_gene_seq);
		this->scenario_n_mismatches+=v_mismatch_list.size();
	}
	if(mismatches_lists.exist(D_gene_seq)){
		const vector<int>& d_mismatch_list = *mismatches_lists[D_gene_seq];
		this->scenario_n_mismatches+=d_mismatch_list.size();
	}
	if(mismatches_lists.exist(J_gene_seq)){
		const vector<int>& j_mismatch_list = *mismatches_lists.at(J_gene_seq);
		this->scenario_n_mismatches+=j_mismatch_list.size();
	}

	if(this->output_scenarios){
		if(this->best_scenarios_vec.size()<this->n_scenarios_counted){

			if(this->best_scenarios_vec.empty()){
				this->best_scenarios_vec.emplace_back(scenario_seq_joint_proba,this->scenario_n_genomic,this->scenario_n_mismatches);
			}
			else{
				vector<tuple<double,size_t,size_t>>::iterator jter = this->best_scenarios_vec.begin();
				while( (scenario_seq_joint_proba>get<0>(*jter)) and (jter != this->best_scenarios_vec.end()) ){
					++jter;
				}
				this->best_scenarios_vec.emplace(jter , scenario_seq_joint_proba,this->scenario_n_genomic,this->scenario_n_mismatches);
			}
		}
		else{
			if(scenario_seq_joint_proba > get<0>(this->best_scenarios_vec[0])){

				auto jter = this->best_scenarios_vec.begin()+1;
				while( (scenario_seq_joint_proba>get<0>(*jter)) and (jter != this->best_scenarios_vec.end()) ){
					++jter;
				}
				this->best_scenarios_vec.emplace(jter , scenario_seq_joint_proba,this->scenario_n_genomic,this->scenario_n_mismatches);
				this->best_scenarios_vec.erase(this->best_scenarios_vec.begin());
			}
		}
	}

	//Add the contributions to the average mutation/genomic counters
	this->sequence_average_n_genomic+=scenario_seq_joint_proba*this->scenario_n_genomic;
	this->sequence_average_n_mismatches+=scenario_seq_joint_proba*this->scenario_n_mismatches;
	this->sequence_average_error_freq+=scenario_seq_joint_proba*((double)this->scenario_n_mismatches/(double)this->scenario_n_genomic);

	//Reset dummy variables
	this->scenario_n_genomic = 0;
	this->scenario_n_mismatches = 0;
}

void Errors_counter::count_sequence(double seq_likelihood , const Model_marginals& single_seq_marginals , const Model_Parms& single_seq_model_parms){
	for(vector<tuple<double,size_t,size_t>>::iterator iter = this->best_scenarios_vec.begin() ; iter!=this->best_scenarios_vec.end() ; ++iter){
		get<0>(*iter)/=seq_likelihood;
		//If an exception is thrown here there is a problem upstream
	}
	//Re normalize the mutation frequency
	this->sequence_average_n_genomic/=seq_likelihood;
	this->sequence_average_n_mismatches/=seq_likelihood;
	this->sequence_average_error_freq/=seq_likelihood;
}

void Errors_counter::initialize_counter(const Model_Parms& parms , const Model_marginals& marginals){
	if(not fstreams_created){
		//Create the individual scenarios file
		if(this->n_scenarios_counted>0){
			output_scenario_errors_file_ptr = shared_ptr<ofstream>(new ofstream);
			this->output_scenario_errors_file_ptr->open(path_to_file + "scenarios_background_and_errors.csv");
			//Create the header
			(*this->output_scenario_errors_file_ptr.get())<<"seq_index;scenario_rank;scenario_proba_cond_seq;n_genomic_nt;n_mismatches"<<endl;
		}

		//Create the average mutation rate scenario file
		output_sequence_averaged_errors_file_ptr = shared_ptr<ofstream>(new ofstream);
		this->output_sequence_averaged_errors_file_ptr->open(path_to_file + "sequence_mutation_frequency.csv");
		//Create the header
		(*this->output_sequence_averaged_errors_file_ptr.get())<<"seq_index;n_genomic_nt;n_mismatches;average_mutation_frequency"<<endl;


		fstreams_created = true;
	}
}

void Errors_counter::add_checked(shared_ptr<Counter> counter){
	return ;
}


void Errors_counter::dump_sequence_data(int seq_index , int iteration_n){

	//Output individual scenarios stats
	if(this->output_scenarios){
		size_t counter = 1;
		for(vector<tuple<double,size_t,size_t>>::reverse_iterator iter = this->best_scenarios_vec.rbegin() ; iter!=this->best_scenarios_vec.rend() ; ++iter){
			(*this->output_scenario_errors_file_ptr.get())<<seq_index<<";"<<counter<<";"<<get<0>(*iter)<<";"<<get<1>(*iter)<<";"<<get<2>(*iter)<<endl;
			++counter;
		}
	}

	//Output sequence stats
	(*this->output_sequence_averaged_errors_file_ptr.get())<<seq_index<<";"<<this->sequence_average_n_genomic<<";"<<this->sequence_average_n_mismatches<<";"<<this->sequence_average_error_freq<<endl;

	best_scenarios_vec.clear();
	this->sequence_average_n_genomic = 0;
	this->sequence_average_n_mismatches = 0;
	this->sequence_average_error_freq = 0;
}

shared_ptr<Counter> Errors_counter::copy() const{
	shared_ptr<Errors_counter> counter_copy_ptr (new Errors_counter(this->n_scenarios_counted));
	counter_copy_ptr->fstreams_created = this->fstreams_created;
	if(this->fstreams_created){
		counter_copy_ptr->output_scenario_errors_file_ptr = this->output_scenario_errors_file_ptr;
		counter_copy_ptr->output_sequence_averaged_errors_file_ptr = this->output_sequence_averaged_errors_file_ptr;
	}
	else{
		throw runtime_error("Counters should not be copied before stream initalization in Errors_counter::copy()");
	}
	return counter_copy_ptr;
}


