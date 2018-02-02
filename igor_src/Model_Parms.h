/*
 * Model_Parms.h
 *
 *  Created on: 3 nov. 2014
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
 *
 *      This class define a Model_Parms object. In our framework this will be used
 *      to easily define the inferred model, by using a directed graph structure.
 *      This helps defining conditional marginals in the model.
 */

#ifndef MODEL_PARMS_H_
#define MODEL_PARMS_H_

#include "Rec_Event.h"
#include "Utils.h"
#include "IntStr.h"
#include <list>
#include <unordered_map>
#include <set>
#include <string>
#include <queue>
#include "Errorrate.h"
#include "Insertion.h"
#include "Deletion.h"
#include "Genechoice.h"
#include "Singleerrorrate.h"
#include "Dinuclmarkov.h"
#include "Hypermutationglobalerrorrate.h"
#include "HypermutationfullNmererrorrate.h"
#include <stdexcept>
#include <memory>

//class Rec_Event;

/**
 * \struct Adjacency_list Model_Parms.h
 * \brief IGoR's Bayesian Network adjacency list
 * \author Q.Marcou
 * \version 1.0
 *
 * Contains a list of smart pointers pointing to an event parents and children (i.e adjacent nodes)
 */
struct Adjacency_list{
	std::list<std::shared_ptr<Rec_Event>> children;
	std::list<std::shared_ptr<Rec_Event>> parents;
};

/**
 * \class Model_Parms Model_Parms.h
 * \brief Implements IGoR's Bayesian Network structure.
 * \author Q.Marcou
 * \version 1.0
 *
 * Implements IGoR's Bayesian Network structure through an acyclic directed graph.
 * Together with the recombination model topology it also contains the error model.
 *	This class implements various methods to extract information from the graph structure such as the order in which RecEvents must be processed provided the topological constraints.
 */
class Model_Parms {
public:
	Model_Parms();
	Model_Parms(std::list<std::shared_ptr<Rec_Event>> event_list);
	Model_Parms(const Model_Parms&);
	//Model_Parms(const Model_Parms&);

	virtual ~Model_Parms();

	//bool is_cyclic(); // adapt an algorithm to find cycles in oriented graphs

	std::list <std::shared_ptr<Rec_Event>> get_children(Rec_Event* ) const;
	std::list<std::shared_ptr<Rec_Event>> get_children(std::shared_ptr<Rec_Event>) const;
	std::list<std::shared_ptr<Rec_Event>> get_children(Rec_Event_name) const;

	std::list <std::shared_ptr<Rec_Event>> get_parents(Rec_Event* ) const;
	std::list <std::shared_ptr<Rec_Event>> get_parents(std::shared_ptr<Rec_Event> ) const;
	std::list <std::shared_ptr<Rec_Event>> get_parents(Rec_Event_name ) const;

	std::list <std::shared_ptr<Rec_Event>> get_ancestors(Rec_Event* ) const;
	std::list <std::shared_ptr<Rec_Event>> get_ancestors(std::shared_ptr<Rec_Event> ) const;
	std::list <std::shared_ptr<Rec_Event>> get_ancestors(Rec_Event_name ) const;

	bool add_edge(Rec_Event* ,Rec_Event*);
	bool add_edge(std::shared_ptr<Rec_Event> , std::shared_ptr<Rec_Event>);
	bool add_edge(Rec_Event_name , Rec_Event_name);

	bool remove_edge(Rec_Event*,Rec_Event*);
	bool remove_edge(std::shared_ptr<Rec_Event>,std::shared_ptr<Rec_Event>);
	bool remove_edge(Rec_Event_name,Rec_Event_name);

	void invert_edge(Rec_Event*,Rec_Event*);
	void invert_edge(std::shared_ptr<Rec_Event>,std::shared_ptr<Rec_Event>);
	void invert_edge(Rec_Event_name,Rec_Event_name);

	bool has_edge(Rec_Event*,Rec_Event*) const;
	bool has_edge(std::shared_ptr<Rec_Event>,std::shared_ptr<Rec_Event>) const;
	bool has_edge(Rec_Event_name,Rec_Event_name) const;

	std::list<std::shared_ptr<Rec_Event>> get_roots() const;

	bool add_event(std::shared_ptr<Rec_Event>);
	bool add_event(Rec_Event*);

	std::queue <std::shared_ptr<Rec_Event>> get_model_queue() const;

	std::shared_ptr<Rec_Event> get_event_pointer(const Rec_Event_name&) const; //const Rec_Event*??
	std::shared_ptr<Rec_Event> get_event_pointer(const std::string& , bool by_nickname) const; //const Rec_Event*??

	void update_edge_event_name(Rec_Event_name,Rec_Event_name);


	void write_model_parms(std::string);

	void read_model_parms(std::string);
	void set_fixed_all_events(bool);

	//Accessors
	std::list<std::shared_ptr<Rec_Event>> get_event_list() const {return events;}

	std::unordered_map<Rec_Event_name,Adjacency_list> get_edges() const {return edges;}

	const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> get_events_map() const;
	std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> get_events_map() ;


	void set_error_ratep(Error_rate* Er_r){error_rate = std::shared_ptr<Error_rate>(Er_r,null_delete<Error_rate>());}
	void set_error_ratep(std::shared_ptr<Error_rate> Er_r){error_rate = Er_r;}

	std::shared_ptr<Error_rate> get_err_rate_p(){return error_rate;}




private:
	std::list <std::shared_ptr<Rec_Event>> events;
	std::unordered_map <Rec_Event_name , Adjacency_list > edges;
	std::shared_ptr<Error_rate> error_rate;


};




#endif /* MODEL_PARMS_H_ */
