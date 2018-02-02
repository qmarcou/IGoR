/*
 * Genechoice.h
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

#ifndef GENECHOICE_H_
#define GENECHOICE_H_

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

/**
 * \class Gene_choice Genechoice.h
 * \brief GeneChoice recombination event.
 * \author Q.Marcou
 * \version 1.0
 *
 * Models the gene choice recombination process.
 * The event realizations are explored based on the sequence alignments that were provdided to the inference.
 * Since D gene can be heavily deleted and might not be recognizable by sequence alignments, a special handling of the D gene choice exploring all D positions ranked by their likelihood has been implemented.
 */
class Gene_choice: public Rec_Event {
	friend class Coverage_err_counter;//Grant friendship to access current gene realization and offset
	friend class Hypermutation_global_errorrate;//Grant friendship to access current gene realization and offset
	friend class Hypermutation_full_Nmer_errorrate;//Same

public:
	//Constructors
	Gene_choice();
	Gene_choice(Gene_class);
	Gene_choice(Gene_class ,std::unordered_map<std::string,Event_realization>&);
	Gene_choice(Gene_class,std::vector<std::pair<std::string,std::string>>);
	//Destructor
	virtual ~Gene_choice();
	//Virtual methods overload
	std::shared_ptr<Rec_Event> copy();
	inline void iterate(double& , Downstream_scenario_proba_bound_map& , const std::string& , const Int_Str& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::shared_ptr<Next_event_ptr>& , Marginal_array_p& , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& , std::shared_ptr<Error_rate>& , std::map<size_t,std::shared_ptr<Counter>>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> & , Safety_bool_map& , Mismatch_vectors_map& , double& , double&);
	void add_realization(int);
	bool add_realization(std::string gene_name , std::string gene_sequence);
	void set_genomic_templates(const std::vector<std::pair<std::string,std::string>>&);
	std::queue<int> draw_random_realization(const Marginal_array_p& , std::unordered_map<Rec_Event_name,int>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::unordered_map<Seq_type , std::string>& , std::default_random_engine&)const;
	void write2txt(std::ofstream&);
	void initialize_event( std::unordered_set<Rec_Event_name>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , Downstream_scenario_proba_bound_map& , Seq_type_str_p_map& , Safety_bool_map& , std::shared_ptr<Error_rate> ,Mismatch_vectors_map&,Seq_offsets_map&,Index_map&);
	void add_to_marginals(long double , Marginal_array_p&) const;


	//Proba bound related computation methods
	bool has_effect_on(Seq_type) const;
	void iterate_initialize_Len_proba( Seq_type considered_junction ,  std::map<int,double>& length_best_proba_map ,  std::queue<std::shared_ptr<Rec_Event>>& model_queue , double& scenario_proba , const Marginal_array_p& model_parameters_point , Index_map& base_index_map , Seq_type_str_p_map& constructed_sequences , int& seq_len ) const ;
	void initialize_Len_proba_bound(std::queue<std::shared_ptr<Rec_Event>>& model_queue , const Marginal_array_p& model_parameters_point , Index_map& base_index_map );



private:
	inline double iterate_common(double ,const int&  , int , Index_map& ,const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& ,const Marginal_array_p&);


	//Inference variables
	//Bool checks
	bool vd_check;
	bool vj_check;
	bool dj_check;

	//Offsets checks


	Seq_Offset d_5_min_offset;
	Seq_Offset d_5_max_offset;
	Seq_Offset j_5_min_offset;
	Seq_Offset j_5_max_offset;
	Seq_Offset v_5_off;
	Seq_Offset v_3_off;
	Seq_Offset d_offset;
	Seq_Offset j_offset;
	Seq_Offset v_offset;




	Seq_Offset v_3_min_offset;
	Seq_Offset v_3_max_offset;
	Seq_Offset d_3_off;
	Seq_Offset d_5_off;



	Seq_Offset d_3_min_offset;
	Seq_Offset d_3_max_offset;


	Seq_Offset j_5_off;

	//Suitable D align bool
	bool no_d_align;
	std::vector<int> no_d_mismatches;
	size_t d_size;
	Seq_Offset d_full_3_offset;



	//Declare common variables
	mutable int base_index;
	double new_scenario_proba;
	double new_tmp_err_w_proba;
	double proba_contribution;
	Int_Str gene_seq;
	int new_index;
	const int* alignment_offset_p;
	std::vector<int>::const_iterator mism_iter;
	std::vector<int>::const_reverse_iterator rev_mism_iter;
	size_t endogeneous_mismatches;


	//Constants
		//Memory Layers
		int memory_layer_cs;
		int memory_layer_mismatches;
		int memory_layer_safety_1;
		int memory_layer_safety_2;
		int memory_layer_off_threep;
		int memory_layer_off_fivep;
		int memory_layer_offset_check1;
		int memory_layer_offset_check2;
		int memory_layer_proba_map_seq;
		int memory_layer_proba_map_junction;
		int memory_layer_proba_map_junction_d2;//If V and J have been chosen D will need to update VJ, VD and DJ
		int memory_layer_proba_map_junction_d3;

		//Gene choices
		bool v_chosen;
		bool v_choice_exist;
		bool d_chosen;
		bool d_choice_exist;
		bool j_chosen;
		bool j_choice_exist;

		//Deletion ranges
		int d_5_max_del;
		int d_5_min_del;
		int d_5_real_max_del;
		int j_5_max_del;
		int j_5_min_del;
		int v_3_max_del;
		int v_3_min_del;
		int d_3_max_del;
		int d_3_min_del;

	//Downstream junction length proba bounds
		std::map<int,double> vd_length_best_proba_map;
		std::map<int,double> vj_length_best_proba_map;
		std::map<int,double> dj_length_best_proba_map;

	//No D prunning proba bound map
		std::map<int,std::vector<std::tuple<std::string,int,int,double>>> vj_length_d_position_proba;

		D_position_comparator D_position_tuple;


};



#endif /* GENECHOICE_H_ */
