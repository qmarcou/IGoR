/*
 * Model_marginals.h
 *
 *  Created on: 3 nov. 2014
 *      Author: marcou
 *      This class defines the model marginals entity, the way they are created and handled.
 */

#ifndef MODEL_MARGINALS_H_
#define MODEL_MARGINALS_H_


//includes from standard libraries
#include <list>
#include <forward_list>
#include <queue>
#include <unordered_map>
#include <stack>
#include <random>
#include <chrono>
#include <vector>
#include <memory>

//Includes from project
#include "Model_Parms.h"
#include "Rec_Event.h"




class Model_marginals {
public:
	Model_marginals();
	Model_marginals(const Model_Parms&);
	Model_marginals(const Model_marginals&);
	virtual ~Model_marginals();
	size_t compute_size(const Model_Parms&);
	size_t get_event_size( std::shared_ptr<const Rec_Event> , const Model_Parms&) const;
	std::pair<std::list<std::pair<Rec_Event_name,size_t>>,std::shared_ptr<double>> compute_event_marginal_probability(Rec_Event_name , const Model_Parms& ) const;
	std::pair<std::list<std::pair<Rec_Event_name,size_t>>,std::shared_ptr<double>> compute_event_marginal_probability(Rec_Event_name , const std::set<Rec_Event_name>& , const Model_Parms& ) const;

	Model_marginals& operator=(const Model_marginals&);
	Model_marginals& operator +=(Model_marginals );
	Model_marginals& operator -=(Model_marginals );
	Model_marginals operator +(Model_marginals );
	Model_marginals operator -(Model_marginals );
	void normalize(std::unordered_map<Rec_Event_name,std::list<std::pair<std::shared_ptr<const Rec_Event>,int>>> , std::unordered_map<Rec_Event_name,int> , std::queue<std::shared_ptr<Rec_Event>>);
	void uniform_initialize(const Model_Parms&);
	void null_initialize();
	void random_initialize(const Model_Parms&);
	void flatten(std::shared_ptr<const Rec_Event>,const Model_Parms& );
<<<<<<< 567d2b3958f347767dadf91d2e59e0d632af4ee0
	void set_realization_proba(std::string,std::shared_ptr<const Rec_Event>,double,const Model_Parms&);
=======
	void add_pseudo_counts(double);
>>>>>>> First attempt for Model_marginals::invert_edge
	bool add_to_marginals(double event_proba , std::list<std::shared_ptr<Rec_Event>> , Model_Parms); //FIXME drop this? pass a list of pointers instead, is this method still needed???
	void copy_fixed_events_marginals(const Model_marginals&,const Model_Parms&,const std::unordered_map<Rec_Event_name,int>&);
	std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>> get_offsets_map(const Model_Parms&) const;
	std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>> get_offsets_map(const Model_Parms& , std::queue<std::shared_ptr<Rec_Event>>) const;
	std::unordered_map<Rec_Event_name,std::list<std::pair<std::shared_ptr<const Rec_Event>,int>>> get_inverse_offset_map(const Model_Parms&) const;
	std::unordered_map<Rec_Event_name,std::list<std::pair<std::shared_ptr<const Rec_Event>,int>>> get_inverse_offset_map(const Model_Parms& , std::queue<std::shared_ptr<Rec_Event>>) const;
	std::unordered_map<Rec_Event_name,int> get_index_map(const Model_Parms&) const; //maybe tie this to the model_marginals itself
	std::unordered_map<Rec_Event_name,int> get_index_map(const Model_Parms&, std::queue<std::shared_ptr<Rec_Event>>) const;
	void write2txt(std::string ,const Model_Parms&);
	void txt2marginals(std::string, const Model_Parms&);
	Model_marginals empty_copy ();
	Model_marginals& invert_edge(Rec_Event_name , Rec_Event_name , Model_Parms&);
	size_t get_length()const{return marginal_arr_size;};


	std::string debug_marg_name;

	//get marginals for given parameter
	std::unique_ptr<long double []> marginal_array_smart_p;

private:
	std::pair<std::list<std::pair<Rec_Event_name,size_t>>,std::shared_ptr<double>> compute_event_marginal_probability(Rec_Event_name , const std::set<Rec_Event_name>& , const Model_Parms& ,const std::unordered_map<Rec_Event_name,int>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , const std::unordered_map<Rec_Event_name,std::list<std::pair<std::shared_ptr<const Rec_Event>,int>>>& ) const;
	void iterate_normalize(std::shared_ptr<const Rec_Event>, std::list<std::pair<std::shared_ptr<const Rec_Event>,int>>& , int ,int );
	void write2txt_iteration(const std::list<std::pair<std::shared_ptr<const Rec_Event>,int>>::const_iterator,const std::list<std::pair<std::shared_ptr<const Rec_Event>,int>>::const_iterator,int,std::ofstream&, std::shared_ptr<Rec_Event> , std::list<std::string>&);
	size_t marginal_arr_size;
	Model_marginals(size_t);

};
struct offset_comp {
	 bool operator()(const std::pair<std::shared_ptr<const Rec_Event>,int> pair_1 , const std::pair<std::shared_ptr<const Rec_Event>,int> pair_2 ){
		 return pair_1.second > pair_2.second;
	 }
};

void recurs_array_copy(std::list<std::shared_ptr<Rec_Event>>::const_iterator,size_t,std::shared_ptr<double>,Marginal_array_p);
void swap_events_order(const Rec_Event_name& ,const Rec_Event_name& , std::pair<std::list<std::pair<Rec_Event_name,size_t>>,std::shared_ptr<double>>&);
void swap_neighboring_events_order(const Rec_Event_name& ,const Rec_Event_name& , std::pair<std::list<std::pair<Rec_Event_name,size_t>>,std::shared_ptr<double>>&);
void align_marginal_array(const std::list<std::pair<Rec_Event_name,size_t>>& , std::pair<std::list<std::pair<Rec_Event_name,size_t>>,std::shared_ptr<double>>&);



#endif /* MODEL_MARGINALS_H_ */
