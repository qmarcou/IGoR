/*
 * Hypermutationglobalerrorrate.cpp
 *
 *  Created on: Mar 8, 2016
 *      Author: quentin
 */

#include "Hypermutationglobalerrorrate.h"

using namespace std;

Hypermutation_global_errorrate::Hypermutation_global_errorrate() {
	// TODO Auto-generated constructor stub

}

Hypermutation_global_errorrate::~Hypermutation_global_errorrate() {
	// TODO Auto-generated destructor stub
}

Error_rate* Hypermutation_global_errorrate::copy()const{

}

Error_rate* Hypermutation_global_errorrate::add_checked(Error_rate* err_r){

}

double Hypermutation_global_errorrate::get_err_rate_upper_bound() const{

}

double Hypermutation_global_errorrate::compare_sequences_error_prob (double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& events_map , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){

}

queue<int> Hypermutation_global_errorrate::generate_errors(string& generated_seq , default_random_engine& generator) const{

}

void Hypermutation_global_errorrate::update(){

	//Update the error rate by gradient descent

}

void Hypermutation_global_errorrate::initialize(const unordered_map<tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& events_map){
	//Get the right pointers for the V gene
	if(learn_on == V_gene | learn_on == VJ_genes | learn_on == VD_genes | learn_on == VDJ_genes){
		try{
			Gene_choice* v_gene_event_p = dynamic_cast<Gene_choice*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side)));
			vgene_offset_p = v_gene_event_p->alignment_offset_p;
			vgene_real_index_p = v_gene_event_p->current_realization_index;

		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize V gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}

	}

	//Get the right pointers for the D gene
	if(learn_on == D_gene | learn_on == DJ_genes | learn_on == VD_genes | learn_on == VDJ_genes){
		try{
			Gene_choice* d_gene_event_p = dynamic_cast<Gene_choice*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side)));
			dgene_offset_p = d_gene_event_p->alignment_offset_p;
			dgene_real_index_p = d_gene_event_p->current_realization_index;
		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize D gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}

	}

	//Get the right pointers for the J gene
	if(learn_on == J_gene | learn_on == DJ_genes | learn_on == VJ_genes | learn_on == VDJ_genes){
		try{
			Gene_choice* j_gene_event_p = dynamic_cast<Gene_choice*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side)));
			jgene_offset_p = j_gene_event_p->alignment_offset_p;
			jgene_real_index_p = j_gene_event_p->current_realization_index;
		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize J gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}
	}
}

void Hypermutation_global_errorrate::add_to_norm_counter(){

}

void Hypermutation_global_errorrate::clean_seq_counters(){

}

void Hypermutation_global_errorrate::write2txt(ofstream& outfile){

}
