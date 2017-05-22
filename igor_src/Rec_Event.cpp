/*
 * Rec_Event.cpp
 *
 *  Created on: 3 nov. 2014
 *      Author: marcou
 */

#include "Rec_Event.h"
#include "Counter.h"

using namespace std;

//std::ofstream log_file(std::string("/media/quentin/419a9e2c-2635-471b-baa0-58a693d04d87/data/tcr_murugan/one_seq_comp/logs.txt"));

Rec_Event::Rec_Event(Gene_class gene , Seq_side side ): priority(0) , event_class(gene) , event_side(side) , name("Undefined_event_name") ,len_min(INT16_MAX) , len_max(INT16_MIN) , type(Undefined_t), event_index(INT16_MIN) , updated(false),fixed(false) , current_realizations_index_vec(vector<int>()) , scenario_downstream_upper_bound_proba(-1),event_upper_bound_proba(-1),scenario_upper_bound_proba(-1),current_realization_index(nullptr){} //FIXME why does this exist? anyway fix initilization


Rec_Event::Rec_Event(Gene_class gene , Seq_side side , unordered_map<string , Event_realization>& realizations): Rec_Event(gene,side)  {
	this->event_realizations = realizations;
}

Rec_Event::Rec_Event(): Rec_Event( Undefined_gene , Undefined_side ) {}



//TODO see this later
/*
Rec_Event::Rec_Event(list<Event_realization> realizations_list){
	for(list<Event_realization>::const_iterator iter = realizations_list.begin() ; iter!=realizations_list.end() ; iter++){
		this->add_realization((*iter));
	}
}


Rec_Event::Rec_Event(list<Event_realization> realization_list, int new_priority ) : Rec_Event(realization_list) {
	this->priority = new_priority;
}
*/


Rec_Event::~Rec_Event() {
	// TODO Auto-generated destructor stub
}

bool Rec_Event::operator ==(const Rec_Event& other)const {
	if(this->get_type() != other.get_type()) return 0;
	if( this->event_class != other.event_class) return 0;
	if( this->event_side != other.event_side) return 0;
	if( this->priority != other.priority) return 0;
	if( this->event_realizations.size() != other.event_realizations.size()) return 0;
	for(unordered_map< string,Event_realization >::const_iterator iter = this->event_realizations.begin() ; iter != this->event_realizations.end() ; ++iter){
		if(other.event_realizations.count((*iter).first) != 1 ) return 0;
	}
	return 1;
}

void Rec_Event::update_event_name(){
	this->name = string() + this->type +string("_")+ this->event_class + string("_") + this->event_side + string("_prio") + to_string(priority) + string("_size") + to_string(this->size());
}

void Rec_Event::add_realization(const Event_realization& realization){
	this->event_realizations.insert( make_pair ((realization).name,(realization)) );
	this->update_event_name();
}



bool Rec_Event::set_priority(int new_priority){
	this->priority = new_priority;
	this->update_event_name();
	return 1;
}


int Rec_Event::size()const{
	return event_realizations.size();
}

void Rec_Event::set_event_identifier(size_t identifier){
	this->event_index = identifier;
}

int Rec_Event::get_event_identifier() const {
	return event_index;
}


void Rec_Event::iterate_wrap_up(double& scenario_proba , Downstream_scenario_proba_bound_map& downstream_proba_map , const std::string& sequence , const Int_Str& int_sequence , Index_map& index_map , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& offset_map , std::shared_ptr<Next_event_ptr>& next_event_ptr_arr  , Marginal_array_p& updated_marginal_array_p , const Marginal_array_p& model_parameters_point ,const std::unordered_map<Gene_class , std::vector<Alignment_data>>& allowed_realizations , Seq_type_str_p_map& constructed_sequences  , Seq_offsets_map& seq_offsets , std::shared_ptr<Error_rate>& error_rate_p , map<size_t,shared_ptr<Counter>>& counters_list ,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& events_map  , Safety_bool_map& safety_set , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){


/*			if(seq_offsets.count(make_pair(J_gene_seq, Three_prime))!=0){

					int offset3=seq_offsets.at(make_pair(J_gene_seq, Three_prime));
					int offset5=seq_offsets.at(make_pair(J_gene_seq, Five_prime));
					if(offset5==0){
						cout<<"problem"<<endl;
					}
					if(offset3==0){
						cout<<"big_problem"<<endl;
					}

	 		  }*/


	if(next_event_ptr_arr.get()[this->event_index]){ //Tests whether the next event pointer is null
			//Recursive call to iterate
			//TODO consider adding a threshold for too low probability events(if necessary)
		next_event_ptr_arr.get()[this->event_index]->iterate(scenario_proba , downstream_proba_map , sequence , int_sequence , index_map , offset_map , next_event_ptr_arr , updated_marginal_array_p , model_parameters_point , allowed_realizations , constructed_sequences  , seq_offsets , error_rate_p , counters_list , events_map , safety_set , mismatches_lists , seq_max_prob_scenario , proba_threshold_factor);


	}
	else{

		long double scenario_error_w_proba = error_rate_p->compare_sequences_error_prob( scenario_proba , sequence , constructed_sequences , seq_offsets , events_map , mismatches_lists , seq_max_prob_scenario , proba_threshold_factor);

		//TODO add a monitor of the likelihood at each iteration

		//TODO implement sequence comparison in the error rate class
		//log_file<<scenario_error_w_proba<<endl;
		//Add the full recombination scenario probability to the marginals
		/*for(forward_list<int*>::const_iterator iter = write_index_list.begin() ; iter!=write_index_list.end() ; iter++){
			updated_marginal_array_p[*(*iter)] += scenario_error_w_proba;
			if(((*iter)==2865) & (!constructed_sequences.at(VD_ins_seq).empty())){
				cout<<"error"<<endl;
			}
			if(((*iter)==2896) & (!constructed_sequences.at(DJ_ins_seq).empty())){
				cout<<"error"<<endl;
			}
		}*/

		if(scenario_error_w_proba>=seq_max_prob_scenario*proba_threshold_factor){
			if(scenario_error_w_proba>seq_max_prob_scenario){seq_max_prob_scenario=scenario_error_w_proba;}

			for(std::map<size_t,std::shared_ptr<Counter>>::iterator iter = counters_list.begin() ; iter != counters_list.end() ; ++iter){
				(*iter).second->count_scenario(scenario_error_w_proba ,scenario_proba , sequence , constructed_sequences , seq_offsets , events_map , mismatches_lists );
			}

			for(std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>::const_iterator iter = events_map.begin() ; iter != events_map.end() ; iter++){
				if(!(*iter).second->is_fixed()){
					(*iter).second->add_to_marginals(scenario_error_w_proba , updated_marginal_array_p);
				}
			}
		}


	}
}



void Rec_Event::initialize_event( unordered_set<Rec_Event_name>& processed_events , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , const unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>>& offset_map , Downstream_scenario_proba_bound_map& downstream_proba_map , Seq_type_str_p_map& constructed_sequences , Safety_bool_map& safety_set , shared_ptr<Error_rate> error_rate_p , Mismatch_vectors_map& mismatches_list , Seq_offsets_map& seq_offsets , Index_map& index_map){
	//No action performed on the event by default if the method is not overloaded
	//Need to call Rec_Event::initialize_event() to apply these common actions when the method is overloaded
	current_realizations_index_vec.push_back(-1);

	if(offset_map.count(this->get_name()) != 0){
		const vector<pair<shared_ptr<const Rec_Event>,int>>& offset_vector = offset_map.at(this->get_name());
		for(vector<pair<shared_ptr<const Rec_Event>,int>>::const_iterator iter = offset_vector.begin() ; iter != offset_vector.end() ; ++iter){
			//Request memory layer
			int event_identitfier = (*iter).first->get_event_identifier();
			index_map.request_memory_layer(event_identitfier);
			memory_and_offsets.emplace_front( event_identitfier , index_map.get_current_memory_layer(event_identitfier) , (*iter).second);
		}
	}

	downstream_proba_map.get_all_current_memory_layer(current_downstream_proba_memory_layers);

	processed_events.emplace(this->name);
	return;
}

void Rec_Event::ind_normalize(Marginal_array_p& marginal_array_p , size_t base_index) const{
	long double sum_marginals = 0;
	for(int i =0 ; i != this->size() ; ++i){
		sum_marginals+= marginal_array_p[base_index + i];
	}
	if(sum_marginals!=0){
		for(int i =0 ; i != this->size() ; ++i){
			marginal_array_p[base_index + i] /= sum_marginals;
		}
	}
}


void Rec_Event::set_crude_upper_bound_proba( size_t base_index , size_t event_size , Marginal_array_p& marginal_array_p){

	double max_proba = 0;
	for(size_t i = 0 ; i!= event_size ; ++i){
		if(marginal_array_p[base_index + i] > max_proba){
			max_proba = marginal_array_p[base_index + i];
		}
	}
	this->event_upper_bound_proba = max_proba;
}

void Rec_Event::set_upper_bound_proba(double proba){
	this->event_upper_bound_proba = proba;
}



/*
 * This method initialize the scenario probability upper bound for each event
 * The point is to compute the upper bound probability (given the model) of the scenario for each event
 * This allows to discard scenarios with too low probability at early stages
 */
void Rec_Event::initialize_crude_scenario_proba_bound(double& downstream_proba_bound , forward_list<double*>& updated_proba_list , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map){
	this->scenario_downstream_upper_bound_proba = downstream_proba_bound;
	this->updated_proba_bounds_list = updated_proba_list;
	if(!this->is_updated()){
		downstream_proba_bound*=this->event_upper_bound_proba;
	}
	else{
		throw logic_error("Updated events should overload Rec_event::initialize_scenario_proba_bound()");
	}
}


/*
 * Description??
 */
double* Rec_Event::get_updated_ptr(){
	throw logic_error("Updated events should overload Rec_event::get_updated_ptr()");
}


/*
 * Updates the value of scenario_upper_bound_proba according to the error weighted scenario and the upper bound of downstream scenarios
 */
void Rec_Event::compute_crude_upper_bound_scenario_proba( double& tmp_err_w_proba ) {
	scenario_upper_bound_proba = tmp_err_w_proba * scenario_downstream_upper_bound_proba;
	for (forward_list<double*>::const_iterator iter = updated_proba_bounds_list.begin() ; iter != updated_proba_bounds_list.end() ; ++iter){
		scenario_upper_bound_proba*=(*(*iter));
	}
}


void Rec_Event::iterate_initialize_Len_proba(Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>>& model_queue , double& scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences ) const{
	int seq_len = 0;
	this->iterate_initialize_Len_proba(considered_junction , length_best_proba_map , model_queue , scenario_proba , model_parameters_point , base_index_map , constructed_sequences , seq_len);
}

/*
 * Called when iterating over all possible scenarios during initialization
 * Fills up the length-max_proba_bound for a given junction , and links the call to iterate_initialize_len_proba for two events
 *
 * TODO constructed sequences should not be used but it is useful to compute the dinucl contribution
 */
void Rec_Event::iterate_initialize_Len_proba_wrap_up(Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>> model_queue , double scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences , int seq_len ) const {

	if(not model_queue.empty()){
		std::shared_ptr<Rec_Event> next_event_p = model_queue.front();
		model_queue.pop();
		//TODO fix this and find a way not to loop over all events
		//if(next_event_p->has_effect_on(considered_junction)){
			// Explore realizations of this event
			next_event_p->iterate_initialize_Len_proba(considered_junction , length_best_proba_map , model_queue , scenario_proba , model_parameters_point , base_index_map , constructed_sequences , seq_len);
		//}
		//else{
			// If this event has no effect on the junction skip it using a recursive call
			//next_event_p->iterate_initialize_Len_proba_wrap_up(considered_junction , length_best_proba_map , model_queue , scenario_proba , model_parameters_point , base_index_map , constructed_sequences , seq_len);
		//}
	}
	else{
		// When all events with an effect on the junction have been processed update the length-proba map
		if(length_best_proba_map.count(seq_len)>0){
			if(scenario_proba>length_best_proba_map.at(seq_len)){
				//Keep the best proba for each length
				length_best_proba_map.at(seq_len) = scenario_proba;
			}
		}
		else{
			length_best_proba_map[seq_len] = scenario_proba;
		}
	}

}

