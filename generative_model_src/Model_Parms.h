/*
 * Model_Parms.h
 *
 *  Created on: 3 nov. 2014
 *      Author: marcou
 *
 *      This class define a Model_Parms object. In our framework this will be used
 *      to easily define the inferred model, by using a directed graph structure.
 *      This helps defining conditional marginals in the model.
 */

#ifndef MODEL_PARMS_H_
#define MODEL_PARMS_H_

#include "Rec_Event.h"
#include "Utils.h"
#include <list>
#include <unordered_map>
#include <string>
#include <queue>
#include "Errorrate.h"
#include "Insertion.h"
#include "Deletion.h"
#include "Genechoice.h"
#include "Singleerrorrate.h"
#include "Dinuclmarkov.h"
#include <stdexcept>

//class Rec_Event;


struct Adjacency_list{
	std::list<Rec_Event*> children;
	std::list<Rec_Event*> parents;
};

class Model_Parms {
public:
	Model_Parms();
	Model_Parms(std::list<Rec_Event*> event_list);
	Model_Parms(const Model_Parms&);
	//Model_Parms(const Model_Parms&);
	virtual ~Model_Parms();
	bool is_cyclic(); // adapt an algorithm to find cycles in oriented graphs
	std::list <Rec_Event*> get_children(Rec_Event* ) const;
	std::list <Rec_Event*> get_parents(Rec_Event* ) const;
	bool add_edge(Rec_Event* ,Rec_Event*);
	bool remove_edge(Rec_Event*,Rec_Event*);
	std::list<Rec_Event*> get_roots() const;
	bool add_event(Rec_Event*);
	std::queue <Rec_Event*> get_model_queue() const;
	Rec_Event* get_event_pointer(const Rec_Event_name&) const; //const Rec_Event*??
	void write_model_parms(std::string);
	void read_model_parms(std::string);

	//Accessors
	std::list<Rec_Event*> get_event_list() const {return events;}
	std::unordered_map<Rec_Event_name,Adjacency_list> get_edges() const {return edges;}
	std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, Rec_Event*> get_events_map();

	void set_error_ratep(Error_rate* Er_r){error_rate = Er_r;}
	Error_rate* get_err_rate_p(){return error_rate;}




private:
	std::list <Rec_Event*> events;
	std::unordered_map <Rec_Event_name , Adjacency_list > edges;
	Error_rate* error_rate;



};

#endif /* MODEL_PARMS_H_ */
