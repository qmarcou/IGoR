/*
 * Bestscenarioscounter.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: quentin
 */

#include "Bestscenarioscounter.h"

using namespace std;

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios): Counter(), n_scenarios_counted(n_scenarios)  {
}

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios , bool is_last_iter_only): Best_scenarios_counter(n_scenarios) {
	this->last_iter_only = is_last_iter_only;
}

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios , string path): Counter(path) , n_scenarios_counted(n_scenarios) {
}

Best_scenarios_counter::Best_scenarios_counter(size_t n_scenarios , string path , bool is_last_iter_only) : Best_scenarios_counter(n_scenarios , path) {
	this->last_iter_only = is_last_iter_only;
}

Best_scenarios_counter::Best_scenarios_counter():Best_scenarios_counter(1,true) {
	// TODO Auto-generated constructor stub

}

Best_scenarios_counter::~Best_scenarios_counter() {
	// TODO Auto-generated destructor stub
}

void Best_scenarios_counter::count_scenario(long double scenario_seq_joint_proba , double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists){


	if(this->best_scenarios_vec.size()<this->n_scenarios_counted){
		for(forward_list<shared_ptr<const Rec_Event>>::const_iterator iter = this->event_fw_list.begin() ; iter != this->event_fw_list.end() ; ++iter){
			this->single_scenario_realizations_queue.push((*iter)->get_current_realizations_index_vec());
		}

		// Get mismatches and add them to the mismatch list
		if(mismatches_lists.exist(V_gene_seq)){
			const vector<int>& v_mismatch_list = *mismatches_lists.at(V_gene_seq);
			single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end() , v_mismatch_list.begin() , v_mismatch_list.end());
		}
		if(mismatches_lists.exist(D_gene_seq)){
			const vector<int>& d_mismatch_list = *mismatches_lists[D_gene_seq];
			single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end() , d_mismatch_list.begin() , d_mismatch_list.end());
		}
		if(mismatches_lists.exist(J_gene_seq)){
			const vector<int>& j_mismatch_list = *mismatches_lists.at(J_gene_seq);
			single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end() , j_mismatch_list.begin() , j_mismatch_list.end());
		}


		if(this->best_scenarios_vec.empty()){
			this->best_scenarios_vec.emplace_back(scenario_seq_joint_proba,const_cast<queue<vector<int>>&>(this->single_scenario_realizations_queue),single_scenario_mismatches_list);
		}
		else{
			vector<tuple<double,queue<vector<int>>,list<int>>>::iterator jter = this->best_scenarios_vec.begin();
			while( (scenario_seq_joint_proba>get<0>(*jter)) and (jter != this->best_scenarios_vec.end()) ){
				++jter;
			}
			this->best_scenarios_vec.emplace(jter , scenario_seq_joint_proba,const_cast<queue<vector<int>>&>(this->single_scenario_realizations_queue),single_scenario_mismatches_list);
		}
	}
	else{
		if(scenario_seq_joint_proba > get<0>(this->best_scenarios_vec[0])){

			for(forward_list<shared_ptr<const Rec_Event>>::const_iterator iter = this->event_fw_list.begin() ; iter != this->event_fw_list.end() ; ++iter){
				this->single_scenario_realizations_queue.push((*iter)->get_current_realizations_index_vec());
			}

			// Get mismatches and add them to the mismatch list
			if(mismatches_lists.exist(V_gene_seq)){
				const vector<int>& v_mismatch_list = *mismatches_lists.at(V_gene_seq);
				single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end() , v_mismatch_list.begin() , v_mismatch_list.end());
			}
			if(mismatches_lists.exist(D_gene_seq)){
				const vector<int>& d_mismatch_list = *mismatches_lists[D_gene_seq];
				single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end() , d_mismatch_list.begin() , d_mismatch_list.end());
			}
			if(mismatches_lists.exist(J_gene_seq)){
				const vector<int>& j_mismatch_list = *mismatches_lists.at(J_gene_seq);
				single_scenario_mismatches_list.insert(single_scenario_mismatches_list.end() , j_mismatch_list.begin() , j_mismatch_list.end());
			}

			auto jter = this->best_scenarios_vec.begin()+1;
			while( (scenario_seq_joint_proba>get<0>(*jter)) and (jter != this->best_scenarios_vec.end()) ){
				++jter;
			}
			this->best_scenarios_vec.emplace(jter , scenario_seq_joint_proba,this->single_scenario_realizations_queue , single_scenario_mismatches_list);
			this->best_scenarios_vec.erase(this->best_scenarios_vec.begin());
		}
	}

	while(not this->single_scenario_realizations_queue.empty()){
		this->single_scenario_realizations_queue.pop();
	}
	single_scenario_mismatches_list.clear();

}

void Best_scenarios_counter::count_sequence(double seq_likelihood , const Model_marginals& single_seq_marginals , const Model_Parms& single_seq_model_parms){
	for(vector<tuple<double,queue<vector<int>>,list<int>>>::iterator iter = this->best_scenarios_vec.begin() ; iter!=this->best_scenarios_vec.end() ; ++iter){
		get<0>(*iter)/=seq_likelihood;
		//If an exception is thrown here there is a problem upstream
	}
}

void Best_scenarios_counter::initialize_counter(const Model_Parms& parms , const Model_marginals& marginals){
	if(not fstreams_created){
		output_scenario_file_ptr = shared_ptr<ofstream>(new ofstream);
		this->output_scenario_file_ptr->open(path_to_file + "best_scenarios_counts.csv");
		//Create the header
		(*this->output_scenario_file_ptr.get())<<"seq_index;scenario_rank;scenario_proba_cond_seq";
		auto event_queue = parms.get_model_queue();

		while(not event_queue.empty()){
			shared_ptr<Rec_Event> event_ptr = event_queue.front();
			(*this->output_scenario_file_ptr.get())<<";"<<event_ptr->get_name();
			this->event_fw_list.emplace_front( event_ptr);
			event_queue.pop();
		}
		event_fw_list.reverse();

		(*this->output_scenario_file_ptr.get())<<";Mismatches"<<endl;

		fstreams_created = true;
	}
	else{
		//Still need to fill in the fw list
		auto event_queue = parms.get_model_queue();
		while(not event_queue.empty()){
			shared_ptr<Rec_Event> event_ptr = event_queue.front();
			this->event_fw_list.emplace_front( event_ptr);
			event_queue.pop();
		}
		event_fw_list.reverse();
	}
}

void Best_scenarios_counter::add_checked(shared_ptr<Counter> counter){
	return ;
}


void Best_scenarios_counter::dump_sequence_data(int seq_index , int iteration_n){

	size_t counter = 1;
	for(vector<tuple<double,queue<vector<int>>,list<int>>>::reverse_iterator iter = this->best_scenarios_vec.rbegin() ; iter!=this->best_scenarios_vec.rend() ; ++iter){
		(*this->output_scenario_file_ptr.get())<<seq_index<<";"<<counter<<";"<<get<0>(*iter);
		queue<vector<int>>& scenario_queue = get<1>(*iter);
		//Loop over events
		while(not scenario_queue.empty()){
			const vector<int>& real_vec = scenario_queue.front();
			(*this->output_scenario_file_ptr.get())<<";(";
			//Loop over event realizations
			for(vector<int>::const_iterator jter = real_vec.begin() ; jter!= real_vec.end() ; ++jter){
				(*this->output_scenario_file_ptr.get())<<(*jter);
				if(jter!=real_vec.end() -1){
					(*this->output_scenario_file_ptr.get())<<",";
				}
			}
			(*this->output_scenario_file_ptr.get())<<")";
			scenario_queue.pop();
		}
		(*this->output_scenario_file_ptr.get())<<";(";
		//Loop over mismatches
		list<int>& mismatches_list = get<2>(*iter);
		list<int>::const_iterator util_iter = mismatches_list.end();
		--util_iter;
		for(list<int>::const_iterator kter = mismatches_list.begin() ; kter!= mismatches_list.end() ; ++kter){
			(*this->output_scenario_file_ptr.get())<<(*kter);
			if(kter!=util_iter){
				(*this->output_scenario_file_ptr.get())<<",";
			}
		}
		(*this->output_scenario_file_ptr.get())<<")"<<endl;

		++counter;
	}

	best_scenarios_vec.clear();
}

shared_ptr<Counter> Best_scenarios_counter::copy() const{
	shared_ptr<Best_scenarios_counter> counter_copy_ptr (new Best_scenarios_counter(this->n_scenarios_counted));
	counter_copy_ptr->fstreams_created = this->fstreams_created;
	if(this->fstreams_created){
		counter_copy_ptr->output_scenario_file_ptr = this->output_scenario_file_ptr;
	}
	else{
		throw runtime_error("Counters should not be copied before stream initalization");
	}
	return counter_copy_ptr;
}

