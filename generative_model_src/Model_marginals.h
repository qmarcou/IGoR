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

//Includes from project
#include "Model_Parms.h"
#include "Rec_Event.h"




class Model_marginals {
public:
	Model_marginals(const Model_Parms&);
	Model_marginals(const Model_marginals&);
	virtual ~Model_marginals();
	size_t compute_size(const Model_Parms&);
	size_t get_event_size( const Rec_Event* , const Model_Parms&) const;
	Model_marginals& operator=(const Model_marginals&);
	Model_marginals& operator +=(Model_marginals );
	Model_marginals operator +(Model_marginals );
	void normalize(std::unordered_map<Rec_Event_name,std::list<std::pair<const Rec_Event*,int>>> , std::unordered_map<Rec_Event_name,int> , std::queue<Rec_Event*>);
	void uniform_initialize(const Model_Parms&);
	void null_initialize();
	void random_initialize(const Model_Parms&);
	void flatten(const Rec_Event*,const Model_Parms& );
	bool add_to_marginals(double event_proba , std::list<Rec_Event*> , Model_Parms); //FIXME drop this? pass a list of pointers instead
	void copy_fixed_events_marginals(const Model_marginals&,const Model_Parms&,const std::unordered_map<Rec_Event_name,int>&);
	std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>> get_offsets_map(const Model_Parms&);
	std::unordered_map<Rec_Event_name,std::vector<std::pair<const Rec_Event*,int>>> get_offsets_map(const Model_Parms& , std::queue<Rec_Event*>);
	std::unordered_map<Rec_Event_name,std::list<std::pair<const Rec_Event*,int>>> get_inverse_offset_map(const Model_Parms&);
	std::unordered_map<Rec_Event_name,std::list<std::pair<const Rec_Event*,int>>> get_inverse_offset_map(const Model_Parms& , std::queue<Rec_Event*>);
	std::unordered_map<Rec_Event_name,int> get_index_map(const Model_Parms&); //maybe tie this to the model_marginals itself
	std::unordered_map<Rec_Event_name,int> get_index_map(const Model_Parms&, std::queue<Rec_Event*>);
	void write2txt(std::string ,const Model_Parms&);
	void txt2marginals(std::string, const Model_Parms&);
	Model_marginals empty_copy ();
	size_t get_length()const{return marginal_arr_size;};



	//get marginals for given parameter
	long double* marginal_array_p;

private:
	void iterate_normalize(Rec_Event*, std::list<std::pair<const Rec_Event*,int>>& , int ,int );
	void write2txt_iteration(const std::list<std::pair<const Rec_Event*,int>>::const_iterator,const std::list<std::pair<const Rec_Event*,int>>::const_iterator,int,std::ofstream&, Rec_Event* , std::list<std::string>&);
	size_t marginal_arr_size;
	Model_marginals(size_t);

};
struct offset_comp {
	 bool operator()(const std::pair<const Rec_Event*,int> pair_1 , const std::pair<const Rec_Event*,int> pair_2 ){
		 return pair_1.second > pair_2.second;
	 }
};



#endif /* MODEL_MARGINALS_H_ */
