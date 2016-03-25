/*
 * Rec_Event.h
 *
 *  Created on: 3 nov. 2014
 *      Author: marcou
 */

#ifndef REC_EVENT_H_
#define REC_EVENT_H_

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <list>
#include <utility>
#include <forward_list>
#include <queue>
#include "Errorrate.h"
#include <random>
#include"Aligner.h"
#include <fstream>
#include "Utils.h"
#include <stdexcept>
#include <tuple>


//class Model_marginals; //forward declare model marginals to avoid circular inclusion

// value of event: struct: event identifier(name of Vgene), event value(sequence), event index(custom)
struct Event_realization {
	const std::string name;
	const int value_int; //union? template? inheritance and reference? just use a virtual class containing two types of events:str and int
	const std::string value_str;
	const std::string value_str_int;
	int index; //Not defined by the user but at the creation of the event, not quite sure about the mutable

	Event_realization(std::string real_name , int val_int , std::string val_str , std::string val_str_int , int index_val): name(real_name) , value_int(val_int) , value_str(val_str) , value_str_int(val_str_int) , index(index_val){}

};








class Rec_Event {
public:
	Rec_Event();
	Rec_Event(Gene_class, Seq_side ,std::unordered_map<std::string,Event_realization>&);
	virtual ~Rec_Event();
	virtual Rec_Event* copy() = 0;//TODO make it const somehow
	virtual int size()const;
	//TODO get rid of deletion map and chosen gene map
	virtual void iterate(double& , double& , const std::string& , const std::string& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>>& , std::queue<Rec_Event*>& , Marginal_array_p& , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& ,Error_rate*&  , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*> & , Safety_bool_map& , Mismatch_vectors_map& , double& , double&)=0 ;
	bool set_priority(int);

	//Accessors
	const Gene_class get_class() const{return event_class;};
	const Seq_side get_side() const{return event_side;};
	const std::unordered_map<std::string , Event_realization>  get_realizations_map() const{return event_realizations;};
	const int get_priority() const{return priority;};
	const Rec_Event_name get_name() const{return name;};
	const std::string get_nickname() const{return nickname;};
	void set_nickname(std::string name){nickname = name;}
	Event_type get_type() const{return this->type;}
	int get_len_max() const {return this->len_max;};
	int get_len_min() const {return this->len_min;};

	bool operator==(const Rec_Event& ) const;
	void update_event_name();
	virtual std::queue<int> draw_random_realization( const Marginal_array_p , std::unordered_map<Rec_Event_name,int>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>>& , std::unordered_map<Seq_type , std::string>& , std::default_random_engine&)const =0;
	virtual void write2txt(std::ofstream&)=0;
	virtual void ind_normalize(Marginal_array_p,size_t);
	virtual void initialize_event( std::unordered_set<Rec_Event_name>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>>& , Seq_type_str_p_map& , Safety_bool_map&  , Error_rate* , Mismatch_vectors_map& , Seq_offsets_map& , Index_map&);
	virtual void initialize_scenario_proba_bound(double& , std::forward_list<double*>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>&);
	virtual void add_to_marginals(long double , Marginal_array_p) const =0;
	virtual void set_upper_bound_proba(size_t , size_t , Marginal_array_p) ;
	void set_upper_bound_proba(double);
	double get_upper_bound_proba()const{return event_upper_bound_proba;};
	//virtual double get_upper_bound_proba() const;
	void set_event_identifier(size_t);
	int get_event_identifier() const;
	bool is_updated() const {return updated;};
	void fix(bool fix_status) {fixed = fix_status;}
	bool is_fixed() const{return fixed;}
	virtual double* get_updated_ptr();
	void compute_upper_bound_scenario_proba( double& ) ;





protected:
	std::unordered_map < std::string, Event_realization > event_realizations;
	int priority;
	Gene_class event_class;
	Seq_side event_side ;
	Rec_Event_name name; //Construct the name in a smart way so that it is unique
	std::string nickname;
	int len_min;
	int len_max;
	Event_type type;
	int event_index;
	std::forward_list<std::tuple<int,int,int>> memory_and_offsets;//0: event identifier , 1: memory layer , 2: offset
	bool updated;
	bool fixed;
	double event_upper_bound_proba;
	double scenario_downstream_upper_bound_proba;
	double scenario_upper_bound_proba; // Used at runtime to store the upper bound probability of the whole scenario
	std::forward_list<double*> updated_proba_bounds_list;


	int compare_sequences(std::string,std::string);//TODO should probably not be a member functino
	void add_realization(const Event_realization&);
	//inline void iterate_wrap_up(double& , double& , const std::string& , const std::string& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>>& , std::queue<Rec_Event*>  , Marginal_array_p&  , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& ,Error_rate*&,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>,const Rec_Event*>&  , Safety_bool_map& , Mismatch_vectors_map& , double& , double&);
	inline void iterate_wrap_up(double& scenario_proba , double& tmp_err_w_proba , const std::string& sequence , const std::string& int_sequence , Index_map& index_map , const std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>>& offset_map , std::queue<Rec_Event*> model_queue  , Marginal_array_p& updated_marginal_array_p , const Marginal_array_p& model_parameters_point ,const std::unordered_map<Gene_class , std::vector<Alignment_data>>& allowed_realizations , Seq_type_str_p_map& constructed_sequences  , Seq_offsets_map& seq_offsets , Error_rate*& error_rate_p,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& events_map  , Safety_bool_map& safety_set , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){


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
				Rec_Event* next_event_p = model_queue.front();
				model_queue.pop();
				//Recursive call to iterate
				//TODO consider adding a threshold for too low probability events(if necessary)
				next_event_p->iterate(scenario_proba , tmp_err_w_proba , sequence , int_sequence , index_map , offset_map , model_queue , updated_marginal_array_p , model_parameters_point , allowed_realizations , constructed_sequences  , seq_offsets , error_rate_p , events_map , safety_set , mismatches_lists , seq_max_prob_scenario , proba_threshold_factor);


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
				for(std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>::const_iterator iter = events_map.begin() ; iter != events_map.end() ; iter++){
					if(!(*iter).second->is_fixed()){
						(*iter).second->add_to_marginals(scenario_error_w_proba , updated_marginal_array_p);
					}
				}
			}


		}
	}





};


//bool compare_events(const Rec_Event*&, const Rec_Event*&);
struct Event_comparator {
	 bool operator()(const std::shared_ptr<Rec_Event> event_p1 , const std::shared_ptr<Rec_Event> event_p2 ){
		 return event_p1->get_priority() > event_p2->get_priority();
	 }
};




#endif /* REC_EVENT_H_ */
