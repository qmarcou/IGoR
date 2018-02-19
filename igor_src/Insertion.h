/*
 * Insertion.h
 *
 *  Created on: Dec 9, 2014
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
 */

#ifndef INSERTION_H_
#define INSERTION_H_

#include "Rec_Event.h"
#include "Utils.h"
#include <forward_list>
#include <unordered_map>
#include <string>
#include <list>
#include <queue>
#include <utility>
#include "Errorrate.h"
#include <random>
#include <map>

/**
 * \class Insertion Insertion.h
 * \brief Insertion recombination events.
 * \author Q.Marcou
 * \version 1.0
 *
 *  The Insertion RecEvent models the distribution of junctional insertion length during the V(D)J recombination process.
 *
 */
class Insertion: public Rec_Event {
public:
	//Constructors
	Insertion();
	Insertion(Gene_class,std::pair<int,int>);
	Insertion(Gene_class,std::forward_list<int>);
	Insertion(Gene_class );
	Insertion(Gene_class  , std::unordered_map<std::string , Event_realization>&);

	//Destructor
	virtual ~Insertion();

	//virtual methods
	std::shared_ptr<Rec_Event> copy();
	inline void iterate(double& , Downstream_scenario_proba_bound_map& , const std::string& , const Int_Str& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::shared_ptr<Next_event_ptr>& , Marginal_array_p& , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& , std::shared_ptr<Error_rate>& , std::map<size_t,std::shared_ptr<Counter>>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> & , Safety_bool_map& , Mismatch_vectors_map& , double& , double&);
	bool add_realization(int);
	std::queue<int> draw_random_realization( const Marginal_array_p& , std::unordered_map<Rec_Event_name,int>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::unordered_map<Seq_type , std::string>& , std::mt19937_64&)const ;
	void write2txt(std::ofstream&);

	void initialize_event( std::unordered_set<Rec_Event_name>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , Downstream_scenario_proba_bound_map& , Seq_type_str_p_map& , Safety_bool_map& , std::shared_ptr<Error_rate> ,Mismatch_vectors_map&,Seq_offsets_map&,Index_map&);
	void add_to_marginals(long double , Marginal_array_p&) const;
	void set_crude_upper_bound_proba(size_t , size_t , Marginal_array_p&) ;
	void initialize_crude_scenario_proba_bound(double& , std::forward_list<double*>& ,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&);

	//Proba bound related computation methods
	bool has_effect_on(Seq_type) const;
	void iterate_initialize_Len_proba( Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>>& model_queue , double& scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences , int& seq_len ) const ;
	void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>>& model_queue , const Marginal_array_p& model_parameters_point , Index_map& base_index_map );

private:
	inline double iterate_common(double , int , int , Index_map& ,const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& ,const Marginal_array_p&);

	std::map<int,Event_realization> ordered_realization_map;
	std::map<int,double> junction_length_best_proba_map;

	mutable Int_Str inserted_str;
	mutable int base_index;
	double new_scenario_proba;
	double proba_contribution;
	int insertions;
	int new_index;
	std::map<int,double> upper_bound_per_ins; //Contains the probability upper bound for each number of insertions
	int previous_index;

	//Iterate common
	int realization_index;
	std::string insertions_str;

	double* dinuc_updated_bound;

	int memory_layer_proba_map_junction;

	//Pre create pairs to call seq_offsets (otherwise cost of creating a pair at each call)
	//std::pair<Seq_type,Seq_side> d_5_pair = std::make_pair (D_gene_seq,Five_prime);
	//std::pair<Seq_type,Seq_side> v_3_pair = std::make_pair (V_gene_seq,Three_prime);
	//std::pair<Seq_type,Seq_side> j_5_pair = std::make_pair (J_gene_seq,Five_prime);
	//std::pair<Seq_type,Seq_side> d_3_pair = std::make_pair (D_gene_seq,Three_prime);




};



#endif /* INSERTION_H_ */
