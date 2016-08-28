/*
 * Insertion.h
 *
 *  Created on: Dec 9, 2014
 *      Author: quentin
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


class Insertion: public Rec_Event {
public:
	//Constructors
	Insertion();
	Insertion(Gene_class,std::pair<int,int>);
	Insertion(Gene_class,std::forward_list<int>);
	Insertion(Gene_class  , std::unordered_map<std::string , Event_realization>&);

	//Destructor
	virtual ~Insertion();

	//virtual methods
	std::shared_ptr<Rec_Event> copy();
	inline void iterate(double& , double& , const std::string& , const Int_Str& , Index_map& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::queue<std::shared_ptr<Rec_Event>>& , Marginal_array_p& , const Marginal_array_p& , const std::unordered_map<Gene_class , std::vector<Alignment_data>>& , Seq_type_str_p_map& , Seq_offsets_map& , std::shared_ptr<Error_rate>& , std::map<size_t,std::shared_ptr<Counter>>& , const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>> & , Safety_bool_map& , Mismatch_vectors_map& , double& , double&);
	bool add_realization(int);
	std::queue<int> draw_random_realization( const Marginal_array_p , std::unordered_map<Rec_Event_name,int>& , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::unordered_map<Seq_type , std::string>& , std::default_random_engine&)const ;
	void write2txt(std::ofstream&);
	void add_to_marginals(long double , Marginal_array_p) const;
	void set_upper_bound_proba(size_t , size_t , Marginal_array_p) ;
	void initialize_scenario_proba_bound(double& , std::forward_list<double*>& ,const std::unordered_map<std::tuple<Event_type,Gene_class,Seq_side>, std::shared_ptr<Rec_Event>>&);

private:
	inline double iterate_common(double , int , int , Index_map& ,const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& ,const Marginal_array_p);

	std::map<int,Event_realization> ordered_realization_map;

	Int_Str inserted_str;
	int base_index;
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

	//Pre create pairs to call seq_offsets (otherwise cost of creating a pair at each call)
	//std::pair<Seq_type,Seq_side> d_5_pair = std::make_pair (D_gene_seq,Five_prime);
	//std::pair<Seq_type,Seq_side> v_3_pair = std::make_pair (V_gene_seq,Three_prime);
	//std::pair<Seq_type,Seq_side> j_5_pair = std::make_pair (J_gene_seq,Five_prime);
	//std::pair<Seq_type,Seq_side> d_3_pair = std::make_pair (D_gene_seq,Three_prime);




};



#endif /* INSERTION_H_ */
