/*
 * Hypermutationglobalerrorrate.cpp
 *
 *  Created on: Mar 8, 2016
 *      Author: quentin
 */

#include "Hypermutationglobalerrorrate.h"

using namespace std;

Hypermutation_global_errorrate::Hypermutation_global_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value){
	//Initialize booleans
	if(apply_to == V_gene | apply_to == VJ_genes | apply_to == VD_genes | apply_to == VDJ_genes){
		apply_to_v = true;
	}
	else apply_to_v = false;

	if(apply_to == D_gene | apply_to == DJ_genes | apply_to == VD_genes | apply_to == VDJ_genes){
		apply_to_d = true;
	}
	else apply_to_d = false;

	if(apply_to == J_gene | apply_to == VJ_genes | apply_to == DJ_genes | apply_to == VDJ_genes){
		apply_to_j = true;
	}
	else apply_to_j = false;

	if(learn_on == V_gene | learn_on == VJ_genes | learn_on == VD_genes | learn_on == VDJ_genes){
		learn_on_v = true;
	}
	else learn_on_v = false;

	if(learn_on == D_gene | learn_on == DJ_genes | learn_on == VD_genes | learn_on == VDJ_genes){
		learn_on_d = true;
	}
	else learn_on_d = false;

	if(learn_on == J_gene | learn_on == VJ_genes | learn_on == DJ_genes | learn_on == VDJ_genes){
		learn_on_j = true;
	}
	else learn_on_j = false;
}

Hypermutation_global_errorrate::~Hypermutation_global_errorrate() {
	// TODO Auto-generated destructor stub
	//Make a clean destructor and delete all the double* contained in maps
}

Error_rate* Hypermutation_global_errorrate::copy()const{

}

Error_rate* Hypermutation_global_errorrate::add_checked(Error_rate* err_r){

}

double Hypermutation_global_errorrate::get_err_rate_upper_bound() const{

}

double Hypermutation_global_errorrate::compare_sequences_error_prob (double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& events_map , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){
	//TODO Take into account the order of mutations

	string& v_gene_seq = (*constructed_sequences[V_gene_seq]);
	string& d_gene_seq = (*constructed_sequences[D_gene_seq]);
	string& j_gene_seq = (*constructed_sequences[J_gene_seq]);

	//First compute the contribution of the errors to the sequence likelihood
	//V gene
	if(apply_to_v){

	}


	//D gene
	if(apply_to_d){

	}


	//J gene
	if(apply_to_j){

	}

	//Record genomic nucleotides coverage and errors

	if(learn_on_v){
		//Get the coverage
		//Get the length of the gene and a pointer to the right array to write on


		//Get the corrected number of deletions
		v_3_del_value_corr = max(0,*v_3_del_value_p);

		// Compute the coverage and nucleotide errors
		//for( i = max(0,-(*vgene_offset_p)) ; i !=  ){

		//}

		//Get the mismatches positions on the gene

	}

	if(learn_on_d){

	}

	if(learn_on_j){

	}

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

			//Initialize gene counters
			const unordered_map<string , Event_realization>& v_realizations = v_gene_event_p->get_realizations_map();
			//Create arrays
			v_gene_nucleotide_coverage_p = new pair<size_t,double*>[v_realizations.size()];
			v_gene_per_nucleotide_error_p = new pair<size_t,double*>[v_realizations.size()];
			v_gene_nucleotide_coverage_seq_p = new pair<size_t,double*>[v_realizations.size()];
			v_gene_per_nucleotide_error_seq_p = new pair<size_t,double*>[v_realizations.size()];

			for(unordered_map<string , Event_realization>::const_iterator iter = v_realizations.begin() ; iter != v_realizations.end() ; iter++){

				//Initialize normalized counters
				v_gene_nucleotide_coverage_p[(*iter).second.index] = pair<size_t,double>((*iter).second.value_str_int.size(),new double [(*iter).second.value_str_int.size()]);
				v_gene_per_nucleotide_error_p[(*iter).second.index] = pair<size_t,double>((*iter).second.value_str_int.size(),new double [(*iter).second.value_str_int.size()]);

				//Initialize sequence counters
				v_gene_nucleotide_coverage_seq_p[(*iter).second.index] = pair<size_t,double>((*iter).second.value_str_int.size(),new double [(*iter).second.value_str_int.size()]);
				v_gene_per_nucleotide_error_seq_p[(*iter).second.index] = pair<size_t,double>((*iter).second.value_str_int.size(),new double [(*iter).second.value_str_int.size()]);
			}

		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize V gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}

		//Get deletion value pointer for V 3' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime)) != 0){
			const Deletion* v_3_del_event_p = dynamic_cast<Deletion*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime)));
			v_3_del_value_p = &(v_3_del_event_p->deletion_value);
		}
		else{v_3_del_value_p = &no_del_buffer;}

	}

	//Get the right pointers for the D gene
	if(learn_on == D_gene | learn_on == DJ_genes | learn_on == VD_genes | learn_on == VDJ_genes){
		try{
			Gene_choice* d_gene_event_p = dynamic_cast<Gene_choice*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side)));
			dgene_offset_p = d_gene_event_p->alignment_offset_p;
			dgene_real_index_p = d_gene_event_p->current_realization_index;
			//Initialize gene counters

		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize D gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}

		//Get deletion value pointer for D 5' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime)) != 0){
			const Deletion* d_5_del_event_p = dynamic_cast<Deletion*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime)));
			d_5_del_value_p = &(d_5_del_event_p->deletion_value);
		}
		else{d_5_del_value_p = &no_del_buffer;}

		//Get deletion value pointer for D 3' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime)) != 0){
			const Deletion* d_3_del_event_p = dynamic_cast<Deletion*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime)));
			d_3_del_value_p = &(d_3_del_event_p->deletion_value);
		}
		else{d_3_del_value_p = &no_del_buffer;}

	}

	//Get the right pointers for the J gene
	if(learn_on == J_gene | learn_on == DJ_genes | learn_on == VJ_genes | learn_on == VDJ_genes){
		try{
			Gene_choice* j_gene_event_p = dynamic_cast<Gene_choice*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side)));
			jgene_offset_p = j_gene_event_p->alignment_offset_p;
			jgene_real_index_p = j_gene_event_p->current_realization_index;
			//Initialize gene counters

		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize J gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}
	}

	//Get deletion value pointer for J 5' deletions if it exists
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime)) != 0){
		const Deletion* j_5_del_event_p = dynamic_cast<Deletion*>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime)));
		j_5_del_value_p = &(j_5_del_event_p->deletion_value);
	}
	else{j_5_del_value_p = &no_del_buffer;}

}

void Hypermutation_global_errorrate::add_to_norm_counter(){

}

void Hypermutation_global_errorrate::clean_seq_counters(){

}

void Hypermutation_global_errorrate::write2txt(ofstream& outfile){

}
