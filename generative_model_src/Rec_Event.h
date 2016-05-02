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
#include "IntStr.h"
#include <random>
#include"Aligner.h"
#include <fstream>
#include "Utils.h"
#include <stdexcept>
#include <tuple>
#include <memory>
#include <map>

class Counter;


//class Model_marginals; //forward declare model marginals to avoid circular inclusion


// value of event: struct: event identifier(name of Vgene), event value(sequence), event index(custom)
struct Event_realization {
	const std::string name;
	const int value_int; //union? template? inheritance and reference? just use a virtual class containing two types of events:str and int
	const std::string value_str;
	const Int_Str value_str_int;
	int index; //Not defined by the user but at the creation of the event, not quite sure about the mutable

	Event_realization(std::string real_name , int val_int , std::string val_str , Int_Str val_str_int , int index_val): name(real_name) , value_int(val_int) , value_str(val_str) , value_str_int(val_str_int) , index(index_val){}

};








class Rec_Event {
public:
	Rec_Event();
	Rec_Event(Gene_class, Seq_side );
	Rec_Event(Gene_class, Seq_side ,std::unordered_map<std::string,Event_realization>&);
	virtual ~Rec_Event();
	virtual std::shared_ptr<Rec_Event> copy() = 0;//TODO make it const somehow
	virtual int size()const;
	//TODO get rid of deletion map and chosen gene map
	virtual void iterate(double& , Downstream_scenario_proba_bound_map& , const std::string& , const Int_Str& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::queue<std::shared_ptr<Rec_Event>>& , Marginal_array_p& , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& , std::shared_ptr<Error_rate>& , std::map<size_t,std::shared_ptr<Counter>>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> & , Safety_bool_map& , Mismatch_vectors_map& , double& , double&)=0 ;
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
	virtual std::queue<int> draw_random_realization( const Marginal_array_p , std::unordered_map<Rec_Event_name,int>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::unordered_map<Seq_type , std::string>& , std::default_random_engine&)const =0;
	virtual void write2txt(std::ofstream&)=0;
	virtual void ind_normalize(Marginal_array_p,size_t);
	virtual void initialize_event( std::unordered_set<Rec_Event_name>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , Downstream_scenario_proba_bound_map& , Seq_type_str_p_map& , Safety_bool_map&  , std::shared_ptr<Error_rate> , Mismatch_vectors_map& , Seq_offsets_map& , Index_map&);
	virtual void initialize_crude_scenario_proba_bound(double& , std::forward_list<double*>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&);
	virtual void add_to_marginals(long double , Marginal_array_p) const =0;
	virtual void set_crude_upper_bound_proba(size_t , size_t , Marginal_array_p) ;
	void set_upper_bound_proba(double);
	double get_upper_bound_proba()const{return event_upper_bound_proba;};
	//virtual double get_upper_bound_proba() const;
	void set_event_identifier(size_t);
	int get_event_identifier() const;
	void set_event_marginal_size(size_t ev_size){this->event_marginal_size = ev_size;};
	bool is_updated() const {return updated;};
	void fix(bool fix_status) {fixed = fix_status;}
	bool is_fixed() const{return fixed;}
	void set_viterbi_run(bool viterbi_like){viterbi_run = viterbi_like;}
	virtual double* get_updated_ptr();
	void compute_crude_upper_bound_scenario_proba( double& ) ;
	const std::vector<int>& get_current_realizations_index_vec() const{return current_realizations_index_vec;};

	//Proba bound related computation methods
	virtual bool has_effect_on(Seq_type) const= 0;
	void iterate_initialize_Len_proba_wrap_up( Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>> model_queue , double scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences , int seq_len ) const;
	virtual void iterate_initialize_Len_proba( Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>>& model_queue , double& scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences , int& seq_len ) const = 0;
	void iterate_initialize_Len_proba( Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>>& model_queue , double& scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences  ) const ;
	virtual void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>>& model_queue , const Marginal_array_p& model_parameters_point , Index_map& base_index_map ) =0;



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
	bool viterbi_run;
	bool initialized;
	size_t event_marginal_size;
	bool fixed;
	double event_upper_bound_proba;
	double scenario_downstream_upper_bound_proba;
	double scenario_upper_bound_proba; // Used at runtime to store the upper bound probability of the whole scenario
	std::forward_list<double*> updated_proba_bounds_list;
	std::vector<int> current_realizations_index_vec;
	const int* current_realization_index;
	int current_downstream_proba_memory_layers[6];





	int compare_sequences(std::string,std::string);//TODO should probably not be a member functino
	void add_realization(const Event_realization&);
	//inline void iterate_wrap_up(double& , double& , const std::string& , const std::string& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>>& , std::queue<Rec_Event*>  , Marginal_array_p&  , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& ,std::shared_ptr<Error_rate>&,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>,const Rec_Event*>&  , Safety_bool_map& , Mismatch_vectors_map& , double& , double&);
	void iterate_wrap_up(double& scenario_proba , Downstream_scenario_proba_bound_map& downstream_proba_map , const std::string& sequence , const Int_Str& int_sequence , Index_map& index_map , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& offset_map , std::queue<std::shared_ptr<Rec_Event>> model_queue  , Marginal_array_p& updated_marginal_array_p , const Marginal_array_p& model_parameters_point ,const std::unordered_map<Gene_class , std::vector<Alignment_data>>& allowed_realizations , Seq_type_str_p_map& constructed_sequences  , Seq_offsets_map& seq_offsets , std::shared_ptr<Error_rate>& error_rate_p , std::map<size_t,std::shared_ptr<Counter>>& counters_list,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& events_map  , Safety_bool_map& safety_set , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor);

};

//bool compare_events(const Rec_Event*&, const Rec_Event*&);
struct Event_comparator {
	 bool operator()(std::shared_ptr<const Rec_Event> event_p1 , std::shared_ptr<const Rec_Event> event_p2 ){
		 return event_p1->get_priority() > event_p2->get_priority();
	 }
};



#endif /* REC_EVENT_H_ */
