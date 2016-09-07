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


Rec_Event::Rec_Event(): event_realizations ( *(new unordered_map<string,Event_realization>)) , priority(0) , event_class(Undefined_gene) , event_side(Undefined_side) , name("Undefined_event_name") , nickname("Undefined_nickname") , len_min(INT16_MAX) , len_max(INT16_MIN) , type(Undefined_t) , event_index(INT16_MIN) , updated(false),fixed(false),current_realizations_index_vec(vector<int>()) {} //FIXME nonsense new

Rec_Event::Rec_Event(Gene_class gene , Seq_side side , unordered_map<string , Event_realization>& realizations): event_realizations(realizations) , priority(0) , event_class(gene) , event_side(side) , name("Undefined_event_name") ,len_min(INT16_MAX) , len_max(INT16_MIN) , type(Undefined_t), event_index(INT16_MIN) , updated(false),fixed(false),current_realizations_index_vec(vector<int>()) {}


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


void Rec_Event::iterate_wrap_up(double& scenario_proba , double& tmp_err_w_proba , const std::string& sequence , const Int_Str& int_sequence , Index_map& index_map , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& offset_map , std::queue<std::shared_ptr<Rec_Event>> model_queue  , Marginal_array_p& updated_marginal_array_p , const Marginal_array_p& model_parameters_point ,const std::unordered_map<Gene_class , std::vector<Alignment_data>>& allowed_realizations , Seq_type_str_p_map& constructed_sequences  , Seq_offsets_map& seq_offsets , std::shared_ptr<Error_rate>& error_rate_p , map<size_t,shared_ptr<Counter>>& counters_list ,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& events_map  , Safety_bool_map& safety_set , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){


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


	if(!model_queue.empty()){
			std::shared_ptr<Rec_Event> next_event_p = model_queue.front();
			model_queue.pop();
			//Recursive call to iterate
			//TODO consider adding a threshold for too low probability events(if necessary)
			next_event_p->iterate(scenario_proba , tmp_err_w_proba , sequence , int_sequence , index_map , offset_map , model_queue , updated_marginal_array_p , model_parameters_point , allowed_realizations , constructed_sequences  , seq_offsets , error_rate_p , counters_list , events_map , safety_set , mismatches_lists , seq_max_prob_scenario , proba_threshold_factor);


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

		for(std::map<size_t,std::shared_ptr<Counter>>::iterator iter = counters_list.begin() ; iter != counters_list.end() ; ++iter){
			(*iter).second->count_scenario(scenario_error_w_proba ,scenario_proba , sequence , constructed_sequences , seq_offsets , events_map , mismatches_lists );
		}


		if(scenario_error_w_proba>=seq_max_prob_scenario*proba_threshold_factor){
			if(scenario_error_w_proba>seq_max_prob_scenario){seq_max_prob_scenario=scenario_error_w_proba;}
			for(std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>::const_iterator iter = events_map.begin() ; iter != events_map.end() ; iter++){
				if(!(*iter).second->is_fixed()){
					(*iter).second->add_to_marginals(scenario_error_w_proba , updated_marginal_array_p);
				}
			}
		}


	}
}



void Rec_Event::initialize_event( unordered_set<Rec_Event_name>& processed_events , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , const unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>>& offset_map , Seq_type_str_p_map& constructed_sequences , Safety_bool_map& safety_set , shared_ptr<Error_rate> error_rate_p , Mismatch_vectors_map& mismatches_list , Seq_offsets_map& seq_offsets , Index_map& index_map){
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

	processed_events.emplace(this->name);
	return;
}

void Rec_Event::ind_normalize(Marginal_array_p marginal_array_p , size_t base_index){
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

void Rec_Event::set_upper_bound_proba( size_t base_index , size_t event_size , Marginal_array_p marginal_array_p){
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
void Rec_Event::initialize_scenario_proba_bound(double& downstream_proba_bound , forward_list<double*>& updated_proba_list , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map){
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
void Rec_Event::compute_upper_bound_scenario_proba( double& tmp_err_w_proba ) {
	scenario_upper_bound_proba = tmp_err_w_proba * scenario_downstream_upper_bound_proba;
	for (forward_list<double*>::const_iterator iter = updated_proba_bounds_list.begin() ; iter != updated_proba_bounds_list.end() ; ++iter){
		scenario_upper_bound_proba*=(*(*iter));
	}
}

