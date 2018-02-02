/*
 * Deletion.h
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
 */

#ifndef DELETION_H_
#define DELETION_H_

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
#include <math.h>

/**
 * \class Deletion Deletion.h
 * \brief Deletion recombination event
 * \author Q.Marcou
 * \version 1.0
 *
 * The Deletion RecEvent models deletions of one genomic fragment on a given side.
 * Deletions can be either positive or negative (= Palindromic insertions)
 *
 * By construction the corresponding GeneChoice must have been explored first.
 */
class Deletion: public Rec_Event {
	friend class Coverage_err_counter;//Grant friendship to access the current number of deletion
	friend class Hypermutation_global_errorrate;//Grant friendship to access the current number of deletion
	friend class Hypermutation_full_Nmer_errorrate;//Same

public:
	//Constructor
	Deletion();
	Deletion(Gene_class, Seq_side ,std::pair<int,int>);
	Deletion(std::forward_list<int>);
	Deletion(Gene_class , Seq_side );
	Deletion(Gene_class , Seq_side , std::unordered_map<std::string , Event_realization>&);

	//Destructor
	virtual ~Deletion();

	//Virtual methods
	std::shared_ptr<Rec_Event> copy();

	inline void iterate(double& , Downstream_scenario_proba_bound_map& , const std::string& , const Int_Str& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::shared_ptr<Next_event_ptr>& , Marginal_array_p& , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& , std::shared_ptr<Error_rate>& , std::map<size_t,std::shared_ptr<Counter>>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> & , Safety_bool_map& , Mismatch_vectors_map& , double& , double&);
	void add_realization(int);
	std::queue<int> draw_random_realization( const Marginal_array_p& , std::unordered_map<Rec_Event_name,int>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::unordered_map<Seq_type , std::string>& , std::default_random_engine&)const;
	void write2txt(std::ofstream&);
	void initialize_event( std::unordered_set<Rec_Event_name>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , Downstream_scenario_proba_bound_map& , Seq_type_str_p_map&  , Safety_bool_map& , std::shared_ptr<Error_rate> , Mismatch_vectors_map&,Seq_offsets_map&,Index_map&);
	void add_to_marginals(long double , Marginal_array_p&) const;




	//Proba bound related computation methods
	bool has_effect_on(Seq_type) const;
	void iterate_initialize_Len_proba( Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>>& model_queue , double& scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences , int& seq_len ) const ;
	void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>>& model_queue , const Marginal_array_p& model_parameters_point , Index_map& base_index_map );


private:
	inline void iterate_common( std::forward_list<Event_realization>::const_iterator& , Index_map& ,const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>&  ,const Marginal_array_p& );

	std::forward_list<Event_realization> int_value_and_index;

	//Declare variables for V part
	bool vd_check;
	bool vj_check;
	Seq_Offset v_3_offset;
	Seq_Offset v_3_new_offset;
	Seq_Offset d_5_min_offset;
	Seq_Offset d_5_max_offset;
	Seq_Offset j_5_min_offset;
	Seq_Offset j_5_max_offset;

	//Declare variables for D part
	//Five prime
	//bool vd_check;
	Seq_Offset d_5_offset ;
	Seq_Offset d_5_new_offset;
	Seq_Offset v_3_min_offset;
	Seq_Offset v_3_max_offset;

	//Three Prime
	bool dj_check;
	Seq_Offset d_3_offset;
	Seq_Offset d_3_new_offset;
	//Seq_Offset j_5_min_offset;
	//Seq_Offset j_5_max_offset;

	//Declare variables for J part
	//bool dj_check;
	//bool vj_check;
	Seq_Offset j_5_offset;
	Seq_Offset j_5_new_offset;
	Seq_Offset d_3_min_offset;
	Seq_Offset d_3_max_offset;
	//Seq_Offset v_3_min_offset;
	//Seq_Offset v_3_max_offset;

	//Initialized variables
	int d_5_max_del;
	int d_5_min_del;
	int j_5_max_del;
	int j_5_min_del;
	int v_3_max_del;
	int v_3_min_del;
	int d_3_max_del;
	int d_3_min_del;

	//Gene choices
	bool v_chosen;
	bool d_chosen;
	bool j_chosen;

	//D_del bool
	bool d_del_opposite_side_processed;

	//Common variables
	mutable int base_index;
	double err_rate_upper_bound;
	double new_scenario_proba;
	double new_tmp_err_w_proba;
	double proba_contribution;
	int new_index;
	//Int_Str previous_str;//&
	mutable Int_Str new_str;
	mutable Int_Str tmp_str;
	mutable std::string gen_new_str;
	mutable std::string gen_tmp_str;
	std::vector<int> mismatches_vector;
	std::vector<int>::iterator mis_iter;
	std::vector<int>::reverse_iterator rev_mis_iter;
	size_t endogeneous_mismatches;
	bool end_reached;
	int deletion_value;

	//Pre create pairs to call seq_offsets (otherwise cost of creating a pair at each call)
	//std::pair<Seq_type,Seq_side> d_5_pair = std::make_pair (D_gene_seq,Five_prime);
	//std::pair<Seq_type,Seq_side> v_3_pair = std::make_pair (V_gene_seq,Three_prime);
	//std::pair<Seq_type,Seq_side> j_5_pair = std::make_pair (J_gene_seq,Five_prime);
	//std::pair<Seq_type,Seq_side> d_3_pair = std::make_pair (D_gene_seq,Three_prime);

	//Seq_type_str_p_map constructed_sequences_copy;
	int memory_layer_cs;
	int memory_layer_mismatches;
	int memory_layer_safety_1;
	int memory_layer_safety_2;
	int memory_layer_offset_del;
	int memory_layer_offset_check1;
	int memory_layer_offset_check2;
	int memory_layer_proba_map_seq;
	int memory_layer_proba_map_junction;

	//Iterate common
	int previous_marginal_index;

	//Downstream junction length proba bounds
		std::map<int,double> vd_length_best_proba_map;
		std::map<int,double> vj_length_best_proba_map;
		std::map<int,double> dj_length_best_proba_map;
};

std::string& make_transversions(std::string& , bool);
Int_Str& make_transversions(Int_Str& );

bool del_numb_compare(const Event_realization& , const Event_realization&) ;




#endif /* DELETION_H_ */
