/*
 * Hypermutationglobalerrorrate.cpp
 *
 *  Created on: Mar 8, 2016
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

#include "Hypermutationglobalerrorrate.h"

using namespace std;

//ofstream debug_stream("/tmp/debug_stream.csv");

/*Hypermutation_global_errorrate::Hypermutation_global_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value): Error_rate() , mutation_Nmer_size(nmer_width) , learn_on(learn) , apply_to(apply) , ei_nucleotide_contributions((new double [4*nmer_width])) , R(starting_flat_value) , n_v_real(0) , n_j_real(0) , n_d_real(0) ,
		v_gene_nucleotide_coverage_p(NULL) , v_gene_per_nucleotide_error_p(NULL),d_gene_nucleotide_coverage_p(NULL) , d_gene_per_nucleotide_error_p(NULL),j_gene_nucleotide_coverage_p(NULL) , j_gene_per_nucleotide_error_p(NULL),
		v_gene_nucleotide_coverage_seq_p(NULL) , v_gene_per_nucleotide_error_seq_p(NULL) , d_gene_nucleotide_coverage_seq_p(NULL) , d_gene_per_nucleotide_error_seq_p(NULL) , j_gene_nucleotide_coverage_seq_p(NULL) , j_gene_per_nucleotide_error_seq_p(NULL) ,
		dj_ins(true) , vd_ins(true) , vj_ins(true) , v_gene(true) , d_gene(true) , j_gene(true) ,
		vgene_offset_p(NULL) , dgene_offset_p(NULL) , jgene_offset_p(NULL) ,
		vgene_real_index_p(NULL) , dgene_real_index_p(NULL) , jgene_real_index_p(NULL),
		v_3_del_value_p(NULL) , d_5_del_value_p(NULL) , d_3_del_value_p(NULL) , j_5_del_value_p(NULL),
		i(-1) , j(-1) , v_3_del_value_corr(INT16_MAX) , d_5_del_value_corr(INT16_MAX) , d_3_del_value_corr(INT16_MAX) , j_5_del_value_corr(INT16_MAX) , tmp_cov_p(NULL) , tmp_err_p(NULL) , tmp_corr_len(-1) , tmp_len_util(-1) , scenario_new_proba(-1) ,
		largest_nuc_adress(-1), tmp_int_nt(-1) , Nmer_index(-1){*/

Hypermutation_global_errorrate::Hypermutation_global_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value): Error_rate() , mutation_Nmer_size(nmer_width) , learn_on(learn) , apply_to(apply) , ei_nucleotide_contributions((new double [4*nmer_width])) , mu(starting_flat_value) , n_v_real(0) , n_j_real(0) , n_d_real(0) ,
		v_sequences(NULL),j_sequences(NULL),
		dj_ins(true) , vd_ins(true) , vj_ins(true) , v_gene(true) , d_gene(true) , j_gene(true) ,
		vgene_offset_p(NULL) , dgene_offset_p(NULL) , jgene_offset_p(NULL) ,
		vgene_real_index_p(NULL) , dgene_real_index_p(NULL) , jgene_real_index_p(NULL),
		v_3_del_value_p(NULL) , d_5_del_value_p(NULL) , d_3_del_value_p(NULL) , j_5_del_value_p(NULL),
		i(-1) , j(-1) , v_3_del_value_corr(INT16_MAX) , d_5_del_value_corr(INT16_MAX) , d_3_del_value_corr(INT16_MAX) , j_5_del_value_corr(INT16_MAX) , tmp_cov_p(NULL) , tmp_err_p(NULL) , tmp_corr_len(-1) , tmp_len_util(-1) , scenario_new_proba(-1) ,
		largest_nuc_adress(-1), tmp_int_nt(-1) , Nmer_index(-1),
		output_Nmer_stat_stream(new ofstream){

	build_upper_bound_matrix(1,1);


	if(fmod(nmer_width,2)==0){
		throw runtime_error("Cannot instanciate hypermutation globale error rate with an even size Nmer(need to be symmetric)");
	}

	size_t array_size = pow(4,mutation_Nmer_size);

	//Instantiate flat nucleotide contributions
	for(int ii=0 ; ii!=4*mutation_Nmer_size ; ++ii){
		ei_nucleotide_contributions[ii] = 0;
	}

	Nmer_mutation_proba = new double [array_size];
	one_seq_Nmer_N_SHM = new double [array_size];
	one_seq_Nmer_N_bg = new double [array_size];
	Nmer_N_SHM = new double [array_size];
	Nmer_N_bg = new double [array_size];
	for(int ii=0 ; ii!=array_size ; ++ii){
		Nmer_mutation_proba[ii] = 0;
		one_seq_Nmer_N_SHM[ii] = 0;
		one_seq_Nmer_N_bg[ii] = 0;
		Nmer_N_SHM[ii] = 0;
		Nmer_N_bg[ii] = 0;
	}

/*	Nmer_mutation_proba = new double [array_size];
	Nmer_P_SHM = new double [array_size];
	Nmer_P_BG = new double [array_size];

	for(int ii=0 ; ii!=array_size ; ++ii){
		Nmer_mutation_proba[ii] = 0;
		Nmer_P_SHM[ii] = 0;
		Nmer_P_BG[ii] = 0;
	}*/


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

	//Initialize adressing vector
	for(int ii = (mutation_Nmer_size-1) ; ii != -1 ; --ii){
		adressing_vector.emplace_back(pow(4,ii));
	}


	empty_vec_util = vector<int>();
	vec_ptr_util = NULL;

	output_Nmer_stat = false;
}

Hypermutation_global_errorrate::Hypermutation_global_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value , vector<double> ei_contributions): Hypermutation_global_errorrate(nmer_width , learn , apply , starting_flat_value){
	if(ei_contributions.size()==(4*mutation_Nmer_size)){
		for(i=0 ; i != ei_contributions.size() ; ++i){
			ei_nucleotide_contributions[i] = ei_contributions[i];
		}
	}
	else{
		throw runtime_error("Size of ei contribution vector does not match the expected size in Hypermutation_global_errorrate(size_t,Gene_class,Gene_class,double,std::vector<double>)");
	}

	this->update_Nmers_proba(0,0,1);
}

Hypermutation_global_errorrate::Hypermutation_global_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value,string filename):Hypermutation_global_errorrate(nmer_width , learn , apply , starting_flat_value) {
	this->set_output_Nmer_stream(filename);
}

Hypermutation_global_errorrate::Hypermutation_global_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value, vector<double> ei_contributions,string filename):Hypermutation_global_errorrate(nmer_width , learn , apply , starting_flat_value,ei_contributions) {
	this->set_output_Nmer_stream(filename);
}

Hypermutation_global_errorrate::~Hypermutation_global_errorrate() {
	// TODO Auto-generated destructor stub
	//Make a clean destructor and delete all the double* contained in maps
	//delete [] ei_nucleotide_contributions;
	//delete [] Nmer_mutation_proba;
	//delete [] Nmer_P_SHM;
	//delete [] Nmer_P_BG;

	//Clean
	if(learn_on_v){
		for(i = 0 ; i != n_v_real ; ++i){
			//delete [] v_gene_nucleotide_coverage_p[i].second;
			//delete [] v_gene_nucleotide_coverage_seq_p[i].second;
			//delete [] v_gene_per_nucleotide_error_p[i].second;
			//delete [] v_gene_per_nucleotide_error_seq_p[i].second;
		}
		//delete [] v_gene_nucleotide_coverage_p; //FIXME find out why free() exception
		//delete [] v_gene_nucleotide_coverage_seq_p;
		//delete [] v_gene_per_nucleotide_error_p;
		//delete [] v_gene_per_nucleotide_error_seq_p;
	}

	if(learn_on_d){

	}

	if(learn_on_j){

	}



}

void Hypermutation_global_errorrate::set_output_Nmer_stream(string filename){
	cout<<"Nmer hypermutation model output set to: "<<filename<<endl;
	output_Nmer_stat_stream->open(filename);
	(*output_Nmer_stat_stream)<<"Nmer_index,N_bg,N_mut"<<endl;
	output_Nmer_stat = true;
}

shared_ptr<Error_rate> Hypermutation_global_errorrate::copy()const{

	shared_ptr<Hypermutation_global_errorrate> copy_err_r = shared_ptr<Hypermutation_global_errorrate>( new Hypermutation_global_errorrate(this->mutation_Nmer_size , this->learn_on , this->apply_to , this->mu) );
	copy_err_r->updated = this->updated;
	copy_err_r->output_Nmer_stat = this->output_Nmer_stat;
	copy_err_r->output_Nmer_stat_stream = this->output_Nmer_stat_stream;
	//copy_err_r->R = this->R;
	for(int ii = 0 ; ii != mutation_Nmer_size*4 ; ++ii){
		copy_err_r->ei_nucleotide_contributions[ii] = this->ei_nucleotide_contributions[ii];
	}
	//Make sure the values for the Nmer probas are correct
	copy_err_r->update_Nmers_proba(0,0,1);

	return copy_err_r;

}

Hypermutation_global_errorrate& Hypermutation_global_errorrate::operator +=(Hypermutation_global_errorrate err_r){

	//FIXME sequential ifs throwing more meaningful exception
	if( (this->learn_on == err_r.learn_on)
		& (this->apply_to == err_r.apply_to)
		& (this->mutation_Nmer_size == err_r.mutation_Nmer_size)
		& (this->mu == err_r.mu)){
		this->number_seq+=err_r.number_seq;
		this->model_log_likelihood+=err_r.model_log_likelihood;


		size_t array_size = pow(4,mutation_Nmer_size);
		for(size_t ii=0 ; ii != array_size ; ++ii){
			//cout<<debug_one_seq_Nmer_N_bg[ii]<<',';
			Nmer_N_SHM[ii] += err_r.Nmer_N_SHM[ii];
			Nmer_N_bg[ii]+= err_r.Nmer_N_bg[ii];

		}



		return *this;
	}
	else{
		throw runtime_error("Hypermutation models cannot be added in Hypermutation_global_errorrate::operator +=()");
	}

}

Error_rate* Hypermutation_global_errorrate::add_checked(Error_rate* err_r){
	return &(this->operator +=( *(dynamic_cast<Hypermutation_global_errorrate*>(err_r) ) ));
}

const double& Hypermutation_global_errorrate::get_err_rate_upper_bound(size_t n_errors , size_t n_error_free) {
/*	double max_proba = 0;
	for(i=0 ; i!=pow(4,mutation_Nmer_size);i++){
		if(Nmer_mutation_proba[i]>max_proba){
			max_proba = Nmer_mutation_proba[i];
		}
	}
	cout<<"max_proba: "<<max_proba<<endl;
	return max_proba/3;*/
	//return mu/(3*(1+mu));
	if( n_errors>this->max_err || n_error_free>this->max_noerr){
		//Need to increase the matrix size (anyway the matrix is at very most read_len^2
		Matrix<double> new_bound_mat (max(this->max_err,n_errors + 10) , max(this->max_noerr , n_error_free+10));
		for(size_t i=0 ; i!=new_bound_mat.get_n_rows() ; ++i){
			for(size_t j=0 ; j!=new_bound_mat.get_n_cols() ; ++j){
				if(i<this->max_err and j<this->max_noerr){
					new_bound_mat(i,j) = this->upper_bound_proba_mat(i,j);
				}
				else{
					new_bound_mat(i,j) = pow(this->mu/(3*(1+mu)),i)*pow(1-mu/((1+mu)),j);
				}
			}
		}
		this->upper_bound_proba_mat = new_bound_mat;
		this->max_err = new_bound_mat.get_n_rows()-1; //-1 since 0 errors is at index 0
		this->max_noerr = new_bound_mat.get_n_cols()-1;
	}
	//TODO find something more sophisticated than mu/(1+mu) approximation

	return this->upper_bound_proba_mat(n_errors,n_error_free);
}

void Hypermutation_global_errorrate::build_upper_bound_matrix(size_t m,size_t n){
	Matrix<double> new_bound_mat (m,n);
	for(size_t i=0 ; i!=new_bound_mat.get_n_rows() ; ++i){
		for(size_t j=0 ; j!=new_bound_mat.get_n_cols() ; ++j){
			if(i<this->max_err and j<this->max_noerr){
				new_bound_mat(i,j) = this->upper_bound_proba_mat(i,j);
			}
			else{
				new_bound_mat(i,j) = pow(this->mu/(3*(1+mu)),i)*pow(1-mu/((1+mu)),j);
			}
		}
	}
	this->upper_bound_proba_mat = new_bound_mat;
	this->max_err = new_bound_mat.get_n_rows()-1;
	this->max_noerr = new_bound_mat.get_n_cols()-1;
}

double Hypermutation_global_errorrate::compare_sequences_error_prob (double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){
	//TODO Take into account the order of mutations?
	//TODO reorganize to be more flexible in the model description

	scenario_resulting_sequence.clear();
	if(v_gene){
		scenario_resulting_sequence += (*constructed_sequences[V_gene_seq]);
	}
	if(d_gene){
		if(vd_ins){
			scenario_resulting_sequence+=(*constructed_sequences[VD_ins_seq]);
		}
		scenario_resulting_sequence+=(*constructed_sequences[D_gene_seq]);
		if(dj_ins){
			scenario_resulting_sequence+=(*constructed_sequences[DJ_ins_seq]);
		}
	}
	else{
		if(vj_ins){
			scenario_resulting_sequence+=(*constructed_sequences[VJ_ins_seq]);
		}
	}
	if(j_gene){
		scenario_resulting_sequence+=(*constructed_sequences[J_gene_seq]);
	}



	vector<int>& v_mismatch_list = *mismatches_lists.at(V_gene_seq);

	if(mismatches_lists.exist(D_gene_seq)){
		vec_ptr_util = mismatches_lists[D_gene_seq];//Should not have to check this and default initialization should be sufficient
	}
	else{
		vec_ptr_util = &empty_vec_util;
	}
	vector<int>& d_mismatch_list = *vec_ptr_util;

	vector<int>& j_mismatch_list = *mismatches_lists.at(J_gene_seq);

	scenario_new_proba = scenario_probability;

	//First compute the contribution of the errors to the sequence likelihood

	//Check that the sequence is at least the Nmer size
	tmp_len_util = scenario_resulting_sequence.size();
	if(tmp_len_util>=mutation_Nmer_size){
		current_mismatch = v_mismatch_list.begin();
		//In case v_mismatch list is empty cannot de-reference afterwards
		if(current_mismatch==v_mismatch_list.end()){
			current_mismatch = d_mismatch_list.begin();
		}
		if(current_mismatch==d_mismatch_list.end()){ //no else if in case d_mismatch list is empty
			current_mismatch = j_mismatch_list.begin();
		}


		//TODO Need to get the previous V nucleotides and last J ones

		//Get the adress of the first Nmer(disregarding the error penalty on the first nucleotides)
		Nmer_index = 0;
		while(!current_Nmer.empty()){
			current_Nmer.pop();
		}

		for(i=0 ; i!=mutation_Nmer_size ; ++i){
			tmp_int_nt = scenario_resulting_sequence.at(i);
			current_Nmer.push(tmp_int_nt);
			Nmer_index+=adressing_vector[i]*tmp_int_nt;
		}
		//FIXME maybe should iterate the other way around, what happens for errors/context of first nucleotides?
		while((current_mismatch!=j_mismatch_list.end())
				&& (*current_mismatch)<(mutation_Nmer_size-1)/2){
			++current_mismatch;

			if(current_mismatch==v_mismatch_list.end()){
				current_mismatch = d_mismatch_list.begin();
			}
			if(current_mismatch==d_mismatch_list.end()){ //no else if in case d_mismatch list is empty
				current_mismatch = j_mismatch_list.begin();
			}

			//Takes care of the fact that current_mismatch is never incremented if there's a mutation at 0 for instance
			//this needs a true correct fix
		}


		//Check if there's an error and apply the cost accordingly

		if( (current_mismatch!=j_mismatch_list.end())
			&& ((*current_mismatch)==(mutation_Nmer_size-1)/2) ){
			scenario_new_proba*=(Nmer_mutation_proba[Nmer_index]/3);
			++current_mismatch;

			if(current_mismatch==v_mismatch_list.end()){
				current_mismatch = d_mismatch_list.begin();
			}
			if(current_mismatch==d_mismatch_list.end()){ //no else if in case d_mismatch list is empty
				current_mismatch = j_mismatch_list.begin();
			}
		}
		else{
			scenario_new_proba*=(1-Nmer_mutation_proba[Nmer_index]);
		}



		//Look at all Nmers in the scenario_resulting_sequence by sliding window
		//Removing the contribution of the first and adding the contribution of the new last

		for( i = (mutation_Nmer_size-1)/2 +1 ; i!=scenario_resulting_sequence.size()-(mutation_Nmer_size-1)/2  ; ++i){

			//Remove the previous first nucleotide of the Nmer and it's contribution to the index
			Nmer_index-=current_Nmer.front()*adressing_vector[0];
			current_Nmer.pop();
			//Shift the index
			Nmer_index*=4;
			//Add the contribution of the new nucleotide
			tmp_int_nt = scenario_resulting_sequence.at(i+(mutation_Nmer_size-1)/2);//Assume a symetric Nmer
			Nmer_index+=tmp_int_nt;
			current_Nmer.push(tmp_int_nt);

			//Apply the error cost
			if( (current_mismatch!=j_mismatch_list.end())
					&& ((*current_mismatch)==i)){

				scenario_new_proba*=(Nmer_mutation_proba[Nmer_index]/3);

				++current_mismatch;

				if(current_mismatch==v_mismatch_list.end()){
					current_mismatch = d_mismatch_list.begin();
				}
				if(current_mismatch==d_mismatch_list.end()){ //no else if in case d_mismatch list is empty
					current_mismatch = j_mismatch_list.begin();
				}
			}
			else{
				if(d_gene){
					if( (i<=seq_offsets.at(V_gene_seq,Three_prime)) or ((i>=seq_offsets.at(D_gene_seq,Five_prime)) and ((i<=seq_offsets.at(D_gene_seq,Three_prime)))) or (i>=seq_offsets.at(J_gene_seq,Five_prime))){
						scenario_new_proba*=(1-Nmer_mutation_proba[Nmer_index]);
						//FIXME THIS A SUPER HARD FIX!
					}
				}
				else{
					if( (i<=seq_offsets.at(V_gene_seq,Three_prime)) or (i>=seq_offsets.at(J_gene_seq,Five_prime))){
						scenario_new_proba*=(1-Nmer_mutation_proba[Nmer_index]);
						//FIXME THIS A SUPER HARD FIX!
					}
				}
			}

		}

///////////////////////////////////////////////////////////////////////////////////////////////////////
		//Hard fix test for taking (mutation_Nmer_size-1)/2 last J nucleotides into account //FIXME
		//This is assuming that J is last nucleotide of the sequence (does agree with the rest of the code until now)
		for(i=scenario_resulting_sequence.size()-(mutation_Nmer_size-1)/2;i!=scenario_resulting_sequence.size();++i){
			Nmer_index-=current_Nmer.front()*adressing_vector[0];
			current_Nmer.pop();
			//Shift the index
			Nmer_index*=4;
			//Add the contribution of the new nucleotide
			tmp_int_nt = j_sequences[**jgene_real_index_p].at(i-seq_offsets.at(J_gene_seq,Five_prime)+(*j_5_del_value_p));//Assume a symetric Nmer
			Nmer_index+=tmp_int_nt;
			current_Nmer.push(tmp_int_nt);

			if( (current_mismatch!=j_mismatch_list.end())
					&& ((*current_mismatch)==i)){

				scenario_new_proba*=(Nmer_mutation_proba[Nmer_index]/3);

				++current_mismatch;

				if(current_mismatch==v_mismatch_list.end()){
					current_mismatch = d_mismatch_list.begin();
				}
				if(current_mismatch==d_mismatch_list.end()){ //no else if in case d_mismatch list is empty
					current_mismatch = j_mismatch_list.begin();
				}
			}
			else{
				if((i>=seq_offsets.at(J_gene_seq,Five_prime))){
					scenario_new_proba*=(1-Nmer_mutation_proba[Nmer_index]);
					//FIXME THIS A SUPER HARD FIX!
				}
			}
		}


////////////////////////////////////////////////////////////////////////////////////////////////
			//Hard fix test for taking (mutation_Nmer_size-1)/2 first V nucleotides into account //FIXME
			//This is assuming that V is the first nucleotide of the sequence (does agree with the rest of the code until now)
		if(v_gene){
			Nmer_index = 0;
			while(!current_Nmer.empty()){
				current_Nmer.pop();
			}

			tmp_len_util = -(**vgene_offset_p) - (mutation_Nmer_size-1)/2 ;
			for(i=0 ; i!=mutation_Nmer_size ; ++i){
				if(i<(mutation_Nmer_size-1)/2){
					//Take an unseen nucleotide from the V
					tmp_int_nt = v_sequences[**vgene_real_index_p][i+tmp_len_util];
				}
				else{
					//Take a nucleotide seen on the read
					tmp_int_nt = scenario_resulting_sequence.at(i-(mutation_Nmer_size-1)/2);
				}
				current_Nmer.push(tmp_int_nt);
				Nmer_index+=adressing_vector[i]*tmp_int_nt;
			}

			//If there is an error on the first nucleotide apply the cost accordingly
			current_mismatch = v_mismatch_list.begin();

			if( (current_mismatch!=v_mismatch_list.end())
				&& ((*current_mismatch)== 0) ){
				scenario_new_proba*=(Nmer_mutation_proba[Nmer_index]/3);
				++current_mismatch;
			}
			else{
				scenario_new_proba*=(1-Nmer_mutation_proba[Nmer_index]);
			}

			for(i=1 ; i!=(mutation_Nmer_size-1)/2;++i){
				Nmer_index-=current_Nmer.front()*adressing_vector[0];
				current_Nmer.pop();
				//Shift the index
				Nmer_index*=4;
				//Add the contribution of the new nucleotide
				tmp_int_nt = scenario_resulting_sequence[i+(mutation_Nmer_size-1)/2];
				Nmer_index+=tmp_int_nt;
				current_Nmer.push(tmp_int_nt);

				if( (current_mismatch!=v_mismatch_list.end())
						&& ((*current_mismatch)==i)){

					scenario_new_proba*=(Nmer_mutation_proba[Nmer_index]/3);

					++current_mismatch;
				}
				else{
					if((i<=seq_offsets.at(V_gene_seq,Three_prime))){
						scenario_new_proba*=(1-Nmer_mutation_proba[Nmer_index]);
						//FIXME THIS A SUPER HARD FIX!
					}
				}
			}
		}

	}



	//If viterbi learning clean seq counters in order to count only this new most likely scenario
	if(viterbi_run){
		this->clean_seq_counters();
	}

	//Record the number of times each Nmer is seen and the number of times it is mutated

	if(learn_on_v){


		/*
		 * If at least (mutation_Nmer_size+1)/2 V nucleotides remaining (1 visible to count the error, (mutation_Nmer_size-1)/2 on the 3' side for the context assessment
		 * => this is most likely ALWAYS the case for V since it is ~300bp long
		 * Length of the visible part of V: seq_offsets.at(V_gene_seq,Three_prime) - seq_offsets.at(V_gene_seq,Five_prime) +1
		 * Length of the non visible part of V: V_gene_size - #visible - #deleted
		 * (could also use vgene offset) => this is what is done for now until I find a good reason why V should not be assumed to be seen all the way 5'
		 * There must be at least one V nucleotide visible (this is ensured by the alignments)
		 *
		 * The situation is simpler for V than J since all non visible nucleotides are included at the beginning
		 */
		if(v_sequences[**vgene_real_index_p].size()- *v_3_del_value_p >= (mutation_Nmer_size+1)/2){
			current_mismatch = v_mismatch_list.begin();

			//Empty the Nmer queue
			Nmer_index = 0;
			while(!current_Nmer.empty()){
				current_Nmer.pop();
			}

			is_visible_nt = false;
			tmp_corr_len = -(**vgene_offset_p) - (mutation_Nmer_size-1)/2;
			tmp_len_util = seq_offsets.at(V_gene_seq,Five_prime)-(mutation_Nmer_size-1)/2; //Start using the information of the (N-1)/2 inserted (or D) nucleotides before the J

			//Fill in the first Nmer queue (=surroundings of the first V nucleotide)
			for(i=0 ; i!=mutation_Nmer_size ; ++i){
				if(is_visible_nt){
					tmp_int_nt = scenario_resulting_sequence.at(i-tmp_len_util);
				}
				else{
					tmp_int_nt = v_sequences[**vgene_real_index_p].at(i+tmp_corr_len);
					if(i==(mutation_Nmer_size-1)/2){
						is_visible_nt = true;
						tmp_len_util = (mutation_Nmer_size-1)/2;
					}

				}
				current_Nmer.push(tmp_int_nt);
				Nmer_index+=adressing_vector.at(i)*tmp_int_nt;
			}

			//Check if there is an error on the first nucleotide and record Nmer statistics
			if( (current_mismatch!=v_mismatch_list.end())
				&& ((*current_mismatch)== 0 )){
				one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				++current_mismatch;
			}
			else{
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
			}

			//Now look at all nucleotides
			/*
			 * i stands for the position of the last nucleotide of the window
			 * Need to stop when i== Vgene 3' offset + #Insertions considered
			 * i.e i == v3' offset + (N-1)/2
			 */
			for(i=(mutation_Nmer_size-1)/2 +1 ; i!= (seq_offsets.at(V_gene_seq,Three_prime) + (mutation_Nmer_size-1)/2 +1) ; ++i){
				//Remove the previous first nucleotide of the Nmer and it's contribution to the index
				Nmer_index-=current_Nmer.front()*adressing_vector[0];
				current_Nmer.pop();
				//Shift the index
				Nmer_index*=4;

				//Get the next int nt
				//In this part all nucleotides are visible
				tmp_int_nt = scenario_resulting_sequence.at(i);

				//Add the contribution of the new nucleotide
				Nmer_index+=tmp_int_nt;
				current_Nmer.push(tmp_int_nt);

				//Check if there is an error on the central nucleotide and record Nmer statistics
				if( (current_mismatch!=v_mismatch_list.end())
						&& ((*current_mismatch)==(i-(mutation_Nmer_size-1)/2))){
					one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
					++current_mismatch;
				}
				else{
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				}
			}

		}
		else{
			//Do nothing: there is not enough nucleotides on the right to compute an Nmer
		}


/*	THIS IS AN OLD VERSION COMMENTED OUT, WILL BE CLEANED SOON
		current_mismatch = v_mismatch_list.begin();

		if(seq_offsets.at(V_gene_seq,Three_prime) >= (mutation_Nmer_size-1)){
			//TODO Need to get the previous V nucleotides and last J ones

			//Get the address of the first Nmer(disregarding the error penalty on the first nucleotides)
			Nmer_index = 0;
			while(!current_Nmer.empty()){
				current_Nmer.pop();
			}

			for(i=0 ; i!=mutation_Nmer_size ; ++i){
				tmp_int_nt = scenario_resulting_sequence.at(i);
				current_Nmer.push(tmp_int_nt);
				Nmer_index+=adressing_vector[i]*tmp_int_nt;
			}

			//FIXME maybe should iterate the other way around, what happens for errors/context of first nucleotides?
			while((current_mismatch!=v_mismatch_list.end())
					&& (*current_mismatch)<(mutation_Nmer_size-1)/2){
				++current_mismatch;
				//Takes care of the fact that current_mismatch is never incremented if there's a mutation at 0 for instance
				//this needs a true correct fix
			}
			//FIXME more!!!!

			//FIXME maybe should iterate the other way around, what happens for errors/context of first nucleotides?

			//Check if there's an error and apply the cost accordingly

			if( (current_mismatch!=v_mismatch_list.end())
				&& ((*current_mismatch)==(mutation_Nmer_size-1)/2) ){
				one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				++current_mismatch;
			}
			else{
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
			}


			if(seq_offsets.at(V_gene_seq,Three_prime) >= (mutation_Nmer_size)){
				//Look at all Nmers in the scenario_resulting_sequence by sliding window
				//Removing the contribution of the first and adding the contribution of the new last

				for( i = (mutation_Nmer_size-1)/2 +1 ; i<seq_offsets.at(V_gene_seq,Three_prime)-(mutation_Nmer_size-1)/2 +1  ; ++i){
					//FIXME seq_offsets.at(V_gene_seq,Three_prime)-(mutation_Nmer_size-1)/2 +1 ??
					//Remove the previous first nucleotide of the Nmer and it's contribution to the index
					Nmer_index-=current_Nmer.front()*adressing_vector[0];
					current_Nmer.pop();
					//Shift the index
					Nmer_index*=4;
					//Add the contribution of the new nucleotide
					tmp_int_nt = scenario_resulting_sequence.at(i+(mutation_Nmer_size-1)/2);//Assume a symetric Nmer
					Nmer_index+=tmp_int_nt;
					current_Nmer.push(tmp_int_nt);

					//Apply the error cost
					if( (current_mismatch!=v_mismatch_list.end())
							&& ((*current_mismatch)==i)){
						one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
						one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
						++current_mismatch;
					}
					else{
						one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
					}


				}
			}
		}*/



	}

	if(learn_on_d){


		if((seq_offsets.at(D_gene_seq,Five_prime)<=seq_offsets.at(D_gene_seq,Three_prime))	//Makes sure there is at least one D nucleotide(not fully deleted)
				and (seq_offsets.at(D_gene_seq,Five_prime)-(mutation_Nmer_size-1)/2 >0) 	//Makes sure there are enough nucleotides on the left
				and (seq_offsets.at(D_gene_seq,Three_prime)+(mutation_Nmer_size-1)/2<scenario_resulting_sequence.size())){	//Makes sure there are enough nucleotides on the right
			current_mismatch = d_mismatch_list.begin();

			//Empty the Nmer queue
			Nmer_index = 0;
			while(!current_Nmer.empty()){
				current_Nmer.pop();
			}

			//tmp_corr_len = seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime)+(mutation_Nmer_size-1)/2;
			tmp_len_util = seq_offsets.at(D_gene_seq,Five_prime)-(mutation_Nmer_size-1)/2; //Start using the information of the (N-1)/2 inserted (or D) nucleotides before the J

			//Fill in the first Nmer queue (=surroundings of the first J nucleotide)
			for(i=0 ; i!= mutation_Nmer_size ; ++i){
				//assume there is no error in the rest of the context => read the scenario resulting sequence
				tmp_int_nt = scenario_resulting_sequence.at(i+tmp_len_util);

				current_Nmer.push(tmp_int_nt);
				Nmer_index+=adressing_vector.at(i)*tmp_int_nt;
			}

			//Check if there is an error on the first nucleotide and record Nmer statistics
			if( (current_mismatch!=d_mismatch_list.end())
				&& ((*current_mismatch)==seq_offsets.at(D_gene_seq,Five_prime)) ){
				one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				++current_mismatch;
			}
			else{
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
			}

			//Now look at all nucleotides
			/*
			 * i stands for the position of the last nucleotide of the window
			 * Need to stop when i== dgene 3' offset + #Insertions/J nucs considered
			 * i.e i == d3' offset + (N-1)/2
			 */
			for(i=(seq_offsets.at(D_gene_seq,Five_prime)+(mutation_Nmer_size-1)/2 +1 ); i!= (seq_offsets.at(D_gene_seq,Three_prime) + (mutation_Nmer_size-1)/2 +1) ; ++i){
				//Remove the previous first nucleotide of the Nmer and it's contribution to the index
				Nmer_index-=current_Nmer.front()*adressing_vector[0];
				current_Nmer.pop();
				//Shift the index
				Nmer_index*=4;

				//Get the next int nt
				//For D all nucleotides are visible
				tmp_int_nt = scenario_resulting_sequence.at(i);

				//Add the contribution of the new nucleotide
				Nmer_index+=tmp_int_nt;
				current_Nmer.push(tmp_int_nt);

				//Check if there is an error on the central nucleotide and record Nmer statistics
				if( (current_mismatch!=d_mismatch_list.end())
						&& ((*current_mismatch)==(i-(mutation_Nmer_size-1)/2))){
					one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
					++current_mismatch;
				}
				else{
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				}
			}



		}
		/*
		 * Some remarks:
		 * could use not visible V and J nucleotides if needed
		 * When D is full deleted d5' offset > d3'offset (in principle at least)
		 */

	}

	if(learn_on_j){

		/*
		 * If at least (mutation_Nmer_size+1)/2 J nucleotides remaining (1 visible to count the error, (mutation_Nmer_size-1)/2 on the 3' side for the context assessment
		 * Length of the visible part of J: seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime) +1
		 * Length of the non visible part of J: J_gene_size - #visible - #deleted
		 * There must be at least one J nucleotide visible (this is ensured by the alignments)
		 */

		if(j_sequences[**jgene_real_index_p].size() - *j_5_del_value_p>=(mutation_Nmer_size+1)/2){
			current_mismatch = j_mismatch_list.begin();


			//Empty the Nmer queue
			Nmer_index = 0;
			while(!current_Nmer.empty()){
				current_Nmer.pop();
			}

			is_visible_nt = true;
			tmp_corr_len = seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime)+(mutation_Nmer_size-1)/2;
			tmp_len_util = seq_offsets.at(J_gene_seq,Five_prime)-(mutation_Nmer_size-1)/2; //Start using the information of the (N-1)/2 inserted (or D) nucleotides before the J

			//Fill in the first Nmer queue (=surroundings of the first J nucleotide)
			for(i=0 ; i!= mutation_Nmer_size ; ++i){
				if(is_visible_nt){
					//For visible nucleotides assume there is no error in the rest of the context => read the scenario resulting sequence
					tmp_int_nt = scenario_resulting_sequence.at(i+tmp_len_util);
					if(i==tmp_corr_len){
						is_visible_nt = false; //All 3' nucleotides are not visible
						tmp_corr_len = (mutation_Nmer_size-1)/2 - *j_5_del_value_p ; //Correct offset to read the j sequence => Should read the J sequence at position i - #insertions + #deletions
					}
				}
				else{
					tmp_int_nt = j_sequences[**jgene_real_index_p].at(i-tmp_corr_len);
				}
				current_Nmer.push(tmp_int_nt);
				Nmer_index+=adressing_vector.at(i)*tmp_int_nt;
			}


			//Check if there is an error on the first nucleotide and record Nmer statistics
			if( (current_mismatch!=j_mismatch_list.end())
				&& ((*current_mismatch)==seq_offsets.at(J_gene_seq,Five_prime)) ){
				one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				++current_mismatch;
			}
			else{
				one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
			}

			//Now look at all nucleotides
			/*
			 * i stands for the position of the last nucleotide of the window
			 * Need to stop when i== #insertions considered + #visible J considered + #invisible J considered
			 * i.e i == (N-1)/2 + (J3'_offset - J5'_offset +1) + (N-1)/2 => i!= N + (J3'_offset - J5'_offset +1)
			 */
			for(i=mutation_Nmer_size ; i!= mutation_Nmer_size + (seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime) +1) ; ++i){
				//Remove the previous first nucleotide of the Nmer and it's contribution to the index
				Nmer_index-=current_Nmer.front()*adressing_vector[0];
				current_Nmer.pop();
				//Shift the index
				Nmer_index*=4;

				//Get the next int nt (either on the visible or invisible part)
				if(is_visible_nt){
					//For visible nucleotides assume there is no error in the rest of the context => read the scenario resulting sequence
					tmp_int_nt = scenario_resulting_sequence.at(i+tmp_len_util);
					if(i==tmp_corr_len){
						is_visible_nt = false; //All 3' nucleotides are not visible
						tmp_corr_len = (mutation_Nmer_size-1)/2 - *j_5_del_value_p; //Correct offset to read the j sequence => #insertion_considered - #deletions + #visible_J_nucs
					}
				}
				else{
					tmp_int_nt = j_sequences[**jgene_real_index_p].at(i-tmp_corr_len);
				}

				//Add the contribution of the new nucleotide
				Nmer_index+=tmp_int_nt;
				current_Nmer.push(tmp_int_nt);

				//Check if there is an error on the central nucleotide and record Nmer statistics
				if( (current_mismatch!=j_mismatch_list.end())
						&& ((*current_mismatch)==i+tmp_len_util-(mutation_Nmer_size-1)/2)){
					one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
					++current_mismatch;
				}
				else{
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				}
			}

		}
		else{
			//Do nothing: there is not enough nucleotides on the right to compute an Nmer
		}


/*	THIS IS AN OLD VERSION COMMENTED OUT, WILL BE CLEANED SOON
 * if( (seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime) +1 )>=mutation_Nmer_size ){
			//There are enough J nucleotides on the read to start counting at the first J position

				current_mismatch = j_mismatch_list.begin();

				//TODO Need to get the previous V nucleotides and last J ones

				//Get the address of the first Nmer(disregarding the error penalty on the first nucleotides)
				Nmer_index = 0;
				while(!current_Nmer.empty()){
					current_Nmer.pop();
				}

				for(i=seq_offsets.at(J_gene_seq,Five_prime)-1 ; i!=seq_offsets.at(J_gene_seq,Five_prime)-1 + mutation_Nmer_size ; ++i){
					tmp_int_nt = scenario_resulting_sequence.at(i);
					current_Nmer.push(tmp_int_nt);
					Nmer_index+=adressing_vector.at(i-seq_offsets.at(J_gene_seq,Five_prime)+1)*tmp_int_nt;
				}


				//FIXME maybe should iterate the other way around, what happens for errors/context of first nucleotides?
				while((current_mismatch!=j_mismatch_list.end())
						&& (*current_mismatch)<(mutation_Nmer_size-1)/2){
					++current_mismatch;
					//Takes care of the fact that current_mismatch is never incremented if there's a mutation at 0 for instance
					//this needs a true correct fix
				}
				//FIXME more!!!!


				//FIXME maybe should iterate the other way around, what happens for errors/context of first nucleotides?

				//Check if there's an error and apply the cost accordingly

				if( (current_mismatch!=j_mismatch_list.end())
					&& ((*current_mismatch)==seq_offsets.at(J_gene_seq,Five_prime)) ){
					one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
					++current_mismatch;
				}
				else{
					one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
				}



				//Look at all Nmers in the scenario_resulting_sequence by sliding window
				//Removing the contribution of the first and adding the contribution of the new last
				if((seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime) )>=mutation_Nmer_size){

					for( i = seq_offsets.at(J_gene_seq,Five_prime) +1 ; i<seq_offsets.at(J_gene_seq,Three_prime)-(mutation_Nmer_size-1)/2 +1  ; ++i){
						//FIXME seq_offsets.at(V_gene_seq,Three_prime)-(mutation_Nmer_size-1)/2 +1 ??
						//Remove the previous first nucleotide of the Nmer and it's contribution to the index
						Nmer_index-=current_Nmer.front()*adressing_vector[0];
						current_Nmer.pop();
						//Shift the index
						Nmer_index*=4;
						//Add the contribution of the new nucleotide
						tmp_int_nt = scenario_resulting_sequence.at(i+(mutation_Nmer_size-1)/2);//Assume a symetric Nmer
						Nmer_index+=tmp_int_nt;
						current_Nmer.push(tmp_int_nt);

						//Apply the error cost
						if( (current_mismatch!=j_mismatch_list.end())
								&& ((*current_mismatch)==i)){
							one_seq_Nmer_N_SHM[Nmer_index] += scenario_new_proba;
							one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
							++current_mismatch;
						}
						else{
							one_seq_Nmer_N_bg[Nmer_index] += scenario_new_proba;
						}


					}
				}
		}
		else if(){
			//There is not enough J nucleotides on the read, however using the "unseen" 3' J nucleotides at least one position can be counted
		}
		else{
			//There is not enough nucleotides on the 3' side to be able to count this position
		}*/


	}


	this->seq_likelihood += scenario_new_proba;
	this->seq_probability+=scenario_probability;
	this->seq_mean_error_number +=  (v_mismatch_list.size() + d_mismatch_list.size() + j_mismatch_list.size())*scenario_new_proba;
	++debug_number_scenarios;

	return scenario_new_proba;

}

queue<int> Hypermutation_global_errorrate::generate_errors(string& generated_seq , mt19937_64& generator) const{
	uniform_real_distribution<double> distribution(0.0,1.0);
	double rand_err ;// distribution(generator);
	queue<int> errors_indices;

	double error_proba;

	Int_Str int_generated_seq = nt2int(generated_seq);

	//FIXME take into account hidden nucleotides on the right and left sides

	//Get the adress of the first Nmer(disregarding the error penalty on the first nucleotides)
	Nmer_index = 0;
	current_Nmer = queue<size_t>();
	for(i=0 ; i!=mutation_Nmer_size ; ++i){
		tmp_int_nt = int_generated_seq.at(i);
		current_Nmer.push(tmp_int_nt);
		Nmer_index+=adressing_vector[i]*tmp_int_nt;
	}

	error_proba = Nmer_mutation_proba[Nmer_index];
	rand_err = distribution(generator);

	if(rand_err<error_proba){
		//Introduce an error
		errors_indices.push((mutation_Nmer_size-1)/2);

		introduce_uniform_transversion(generated_seq[(mutation_Nmer_size-1)/2], generator , distribution);
	}

	for( i = (mutation_Nmer_size+1)/2 ; i!=int_generated_seq.size()-(mutation_Nmer_size-1)/2 ; ++i){
		//Remove the previous first nucleotide of the Nmer and it's contribution to the index
		Nmer_index-=current_Nmer.front()*adressing_vector[0];
		current_Nmer.pop();
		//Shift the index
		Nmer_index*=4;
		//Add the contribution of the new nucleotide
		tmp_int_nt = int_generated_seq.at(i+(mutation_Nmer_size-1)/2);//Assume a symmetrically sized Nmer
		Nmer_index+=tmp_int_nt;
		current_Nmer.push(tmp_int_nt);


		error_proba = Nmer_mutation_proba[Nmer_index];
		rand_err = distribution(generator);

		if(rand_err<error_proba){
			//Introduce an error
			errors_indices.push(i);

			introduce_uniform_transversion(generated_seq[i], generator , distribution);
		}
	}
	return errors_indices;
}

uint64_t Hypermutation_global_errorrate::generate_random_contributions(double ei_contribution_range){
	//Create seed for random generator
	//create a seed from timer
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point time = myclock::now();
	myclock::duration dur = myclock::time_point::max() - time;

	//Get a random seed
	uint64_t random_seed = draw_random_64bits_seed();
	//Instantiate random number generator
	mt19937_64 generator =  mt19937_64(random_seed);
	uniform_real_distribution<double> distribution(-ei_contribution_range,ei_contribution_range);

	for(i = 0 ; i != mutation_Nmer_size ; ++i){
		double contribution_sum = 0;
		double rand_contribution;
		for(j = 0 ; j != 4 ; ++j){
			rand_contribution = distribution(generator);
			ei_nucleotide_contributions[i*4+j] = rand_contribution;
			contribution_sum += rand_contribution;
		}
		//Ensure that the sum of the contributions at one position is 0
		for(j=0 ; j != 4 ; ++j){
			cout<<ei_nucleotide_contributions[i*4+j]<<";";
			ei_nucleotide_contributions[i*4+j] -= contribution_sum/4;
			cout<<ei_nucleotide_contributions[i*4+j]<<endl;
		}
	}
	this->update_Nmers_proba(0,0,1);

	return random_seed;
}


void Hypermutation_global_errorrate::update(){

	//If an output stream is provided outputs Nmer statistics in a file (mostly for debugging)
	if(output_Nmer_stat){
		for(int zzz=0 ; zzz!=pow(4,mutation_Nmer_size) ; ++zzz){
			(*output_Nmer_stat_stream)<<zzz<<","<<Nmer_N_bg[zzz]<<","<<Nmer_N_SHM[zzz]<<endl;
		}
	}


	//Update the error rate by maximizing the likelihood of the error model

	/*
	 * Find the maximum likelihood parameters by using newton's method on the derivative
	 * The constraint on the sum of e_j at one position being zero is absorbed by setting
	 * e_j(\pi_k) = -sum(e_j(\pi_i)) for i!=k
	 * At each step we solve H\deltaX=-J (where H is the hessian and J the jacobian of the "modified error-likelihood" function)
	 */

	double j_norm = INFINITY;
	double x_norm = INFINITY;

	size_t counter = 0;

	while(x_norm>1e-5 && counter<50000){

		double error_model_likelihood = 0;

		//Construct the 3N+1 sized Jacobian vector
		double J_data[3*mutation_Nmer_size+1];
		for(i=0 ; i!= 3*mutation_Nmer_size+1 ; ++i){
			J_data[i] = 0;
		}
		gsl_vector_view J = gsl_vector_view_array (J_data, 3*mutation_Nmer_size+1);

		//Construct the 3N+1 square Hessian matrix
		double H_data[(3*mutation_Nmer_size+1)*(3*mutation_Nmer_size+1)]; //Are gsl matrices row or column first?
		for(i=0 ; i!= (3*mutation_Nmer_size+1)*(3*mutation_Nmer_size+1) ; ++i){
			H_data[i] = 0;
		}
		gsl_matrix_view H = gsl_matrix_view_array (H_data, 3*mutation_Nmer_size+1,3*mutation_Nmer_size+1);

/*		for(int yyy = 0 ; yyy !=(3*mutation_Nmer_size)+1 ; ++yyy){
			for(int zzz = 0 ; zzz !=(3*mutation_Nmer_size)+1 ; ++zzz){
				cout<<gsl_matrix_get(&H.matrix , zzz , yyy)<<";";
			}
			cout<<endl;
		}
		cout<<endl;*/

		//Compute the values for the Jacobian and Hessian entries
		int base_4_address[mutation_Nmer_size];
		int max_address = 0;
		for (i=0;i!=mutation_Nmer_size;++i){
			base_4_address[i]=0;
			max_address += 3*adressing_vector[i];
		}
		max_address+=1;
		j=0;
		while(j!=max_address){
			//double current_Nmer_P_SHM = Nmer_P_SHM[j];
			//double current_Nmer_P_bg = Nmer_P_BG[j];
			double current_Nmer_P_SHM = Nmer_N_SHM[j];
			double current_Nmer_P_bg = Nmer_N_bg[j];
			double current_Nmer_unorm_score = compute_Nmer_unorm_score(base_4_address,ei_nucleotide_contributions);
			if(current_Nmer_P_bg!=0){
				for(i=0;i!=mutation_Nmer_size;++i){
					if(base_4_address[i]==3){
						//Add contribution to the 3 constraining nucleotides
						for(size_t jj=0 ; jj!=3 ; ++jj){

							//Add contribution to dQ/dei
							J_data[i*3 + jj] -= current_Nmer_P_SHM - current_Nmer_P_bg*(mu*current_Nmer_unorm_score)/(1+mu*current_Nmer_unorm_score) ;

							//Add contribution to d²Q/dei²
							*(gsl_matrix_ptr( &H.matrix, i*3 + jj,i*3 + jj)) += -current_Nmer_P_bg *(mu*current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);

							//Add contributions to the two other constrained nucleotides
							for(size_t jjj=(jj+1) ; jjj != 3 ; ++jjj){
								*(gsl_matrix_ptr( &H.matrix, i*3 + jj,i*3 + jjj)) += -current_Nmer_P_bg *(mu*current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);
							}

							//Add contribution to d²Q/dejdei
							for(size_t ii=(i+1) ; ii!=mutation_Nmer_size ; ++ii){
								//Since the Hessian is symmetric no need to go over everything twice(only fill the lower triangular matrix)
								if(base_4_address[ii]==3){
									for(size_t jjj=0 ; jjj!=3 ; ++jjj){
										*(gsl_matrix_ptr(&H.matrix,i*3 + jj ,ii*3 + jjj)) += -current_Nmer_P_bg *(mu*current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);
									}
								}
								else{
									*(gsl_matrix_ptr(&H.matrix , i*3 + jj  , ii*3 + base_4_address[ii])) -= -current_Nmer_P_bg *(mu*current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);
								}
							}

							//Add contribution to d²Q/dRdei
							*(gsl_matrix_ptr( &H.matrix , i*3 + jj , 3*mutation_Nmer_size)) -= -current_Nmer_P_bg *(current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);

						}
					}
					else{
						//Add contribution to dQ/dei
						J_data[i*3 + base_4_address[i]] += current_Nmer_P_SHM - current_Nmer_P_bg*(mu*current_Nmer_unorm_score)/(1+mu*current_Nmer_unorm_score) ;

						//Add contribution to d²Q/dei²
						*(gsl_matrix_ptr( &H.matrix ,i*3 + base_4_address[i] , i*3 + base_4_address[i])) += -current_Nmer_P_bg *(mu*current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);


						//Add contribution to d²Q/dejdei
						for(size_t ii=(i+1) ; ii!=mutation_Nmer_size ; ++ii){
							//Since the Hessian is symmetric no need to go over everything twice(only fill the lower triangular matrix)
							if(base_4_address[ii]==3){
								for(size_t jj=0 ; jj!=3 ; ++jj){
									*(gsl_matrix_ptr( &H.matrix , i*3 + base_4_address[i],ii*3 + jj )) -= -current_Nmer_P_bg *(mu*current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);
								}
							}
							else{
								*(gsl_matrix_ptr( &H.matrix , i*3 + base_4_address[i],ii*3 + base_4_address[ii] )) += -current_Nmer_P_bg *(mu*current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);
							}

							//if same site(ii==i) the only remaining term is the one from the constrained contribution(4th nucleotide), taken care of in the upper part
						}
						//Add contribution to d²Q/dRdei
						*(gsl_matrix_ptr( &H.matrix , i*3 + base_4_address[i] , 3*mutation_Nmer_size)) += -current_Nmer_P_bg *(current_Nmer_unorm_score)/pow(1+mu*current_Nmer_unorm_score,2);
					}



				}
				//Add contribution to dQ/dR and d²Q/dR²
				J_data[(3*mutation_Nmer_size)] += current_Nmer_P_SHM/mu -(current_Nmer_P_bg*current_Nmer_unorm_score/(1+mu*current_Nmer_unorm_score));
				H_data[(3*mutation_Nmer_size+1)*(3*mutation_Nmer_size+1)-1] += current_Nmer_P_bg*pow(current_Nmer_unorm_score,2)/(1+mu*current_Nmer_unorm_score) -current_Nmer_P_SHM/pow(mu,2);

				//Copy the symmetric part of the Hessian matrix
				for(int ii = 0 ; ii!=3*mutation_Nmer_size+1 ; ++ii){
					for(int jj = ii+1 ; jj!=3*mutation_Nmer_size+1 ; ++jj){
						*(gsl_matrix_ptr( &H.matrix , jj , ii)) = *(gsl_matrix_ptr( &H.matrix , ii , jj));
					}
				}//FIXME is this supposed to be here?

			}

			//Update the base 10 and 4 addresses
			this->increment_base_10_and_4(j,base_4_address);

			error_model_likelihood+=current_Nmer_P_SHM*(log(mu)+log(current_Nmer_unorm_score) - log(3)) - current_Nmer_P_bg*log(1+mu*current_Nmer_unorm_score);
		}

		//cout<<"current hypermutation model likelihood: "<<error_model_likelihood<<endl;

		j_norm = 0;
		for(int kk=0 ; kk != (3*mutation_Nmer_size+1) ; ++kk){
			j_norm += pow(J_data[kk],2);
		}
		j_norm = sqrt(j_norm);
		if(j_norm==0){
			cout<<"Newton's method converged"<<endl;
			break;
		}
		else if(std::isnan(j_norm)){
			throw runtime_error("Optimization of the hypermutation model failed");
		}
		else{
			cout<<"jacobian norm: "<<j_norm<<endl;
		}

		//Set J to -J and then solve H\deltax = J
		for(i=0 ; i!= 3*mutation_Nmer_size+1 ; ++i){
			J_data[i] = -J_data[i];
		}

/*		cout<<"Jacobian"<<endl;
		for(int zzz = 0 ; zzz !=(3*mutation_Nmer_size)+1 ; ++zzz){
			cout<<gsl_vector_get(&J.vector , zzz )<<";";
		}
		cout<<endl;

		cout<<"Hessian"<<endl;
		for(int yyy = 0 ; yyy !=(3*mutation_Nmer_size)+1 ; ++yyy){
			for(int zzz = 0 ; zzz !=(3*mutation_Nmer_size)+1 ; ++zzz){
				cout<<gsl_matrix_get(&H.matrix , zzz , yyy)<<";";
			}
			cout<<endl;
		}
		cout<<endl;*/

		//Solve the system
		gsl_vector *x = gsl_vector_alloc (3*mutation_Nmer_size+1);
		gsl_permutation * p = gsl_permutation_alloc (3*mutation_Nmer_size+1);
		int signum;
		gsl_linalg_LU_decomp(&H.matrix, p , &signum);
		gsl_linalg_LU_solve (&H.matrix, p, &J.vector, x);

		//solve the system using a pseudo inverse
/*		gsl_matrix *V = gsl_matrix_alloc(3*mutation_Nmer_size+1,3*mutation_Nmer_size+1);
		gsl_vector *S = gsl_vector_alloc(3*mutation_Nmer_size+1);
		gsl_vector *work_vect = gsl_vector_alloc(3*mutation_Nmer_size+1);
		gsl_linalg_SV_decomp(&H.matrix , V , S , work_vect);

		gsl_linalg_SV_solve (&H.matrix, V , S , &J.vector, x );	*/

/*		cout<<"\\deltaX vector"<<endl;
		for(int ii= 0 ; ii != mutation_Nmer_size*3+1 ; ++ii){
			cout<<gsl_vector_get(x,ii)<<";";
		}
		cout<<endl;*/

		//Make x a unit direction vector
		x_norm =  gsl_blas_dnrm2 ( x );
		gsl_vector_scale (x , 1.0/x_norm );

		//Backtracking line search
		double alpha = x_norm;
		double tau=.999;
		double c=.01;

		double m;



		//Get the dot product of this direction and the Gradient
		//Reset J to -J(get the gradient)
		for(i=0 ; i!= 3*mutation_Nmer_size+1 ; ++i){
			J_data[i] = -J_data[i];
		}
		gsl_blas_ddot ( x , &J.vector , &m );

		//cout<<"start backtracking line search"<<endl;
		//cout<<"m:"<<m<<endl;
		//cout<<"x_norm"<<gsl_blas_dnrm2 ( x )<<endl;
		//cout<<"grad_norm: "<<gsl_blas_dnrm2 ( &J.vector )<<endl;
		//cout<<alpha<<endl;
		//Now reduce alpha until Armijo–Goldstein condition is fulfilled
		double new_error_model_likelihood = compute_new_model_likelihood(alpha,x);
		while(new_error_model_likelihood-error_model_likelihood< alpha*c*m ){
			alpha*=tau;
			new_error_model_likelihood = compute_new_model_likelihood(alpha,x);
			//cout<<"alpha: "<<alpha<<", ML:"<<new_error_model_likelihood-error_model_likelihood<<";"<<new_error_model_likelihood<<";"<<error_model_likelihood<<endl;
			//std::cin.ignore();
		}
		//cout<<alpha<<endl;


		//Update the parameters values
		for(i=0 ; i != mutation_Nmer_size ; ++i){
			//ei_nucleotide_contributions[i*mutation_Nmer_size+3] = 0;
			for(j=0 ; j!=3 ; ++j){
				//Update the contribution of the nucleotide
				ei_nucleotide_contributions[i*4+j] += alpha*gsl_vector_get(x,(i*3 + j));

				//Compute the contribution of the constrained nucleotide
				ei_nucleotide_contributions[i*4+3] -= alpha*gsl_vector_get(x,(i*3 + j));
			}
		}
		//Update the normalization factor
		mu += alpha*gsl_vector_get(x,(3*mutation_Nmer_size));

/*		cout<<"new model parms"<<endl;
		cout<<R<<endl;
		for(int zzz=0 ; zzz!=4*mutation_Nmer_size ; ++zzz){
			cout<<ei_nucleotide_contributions[zzz]<<";";
		}
		cout<<endl;
		cout<<endl;*/

		++counter;

	}

	//Compute the new mutation probabilities for the full Nmers
	this->update_Nmers_proba(0,0,1);

	//Clean counters
	this->clean_all_counters();

}

double Hypermutation_global_errorrate::compute_new_model_likelihood(double alpha ,gsl_vector* update_vect_p){
	double new_R = mu + alpha*gsl_vector_get(update_vect_p,(3*mutation_Nmer_size));
	double new_ei_nucleotide_contributions [4*mutation_Nmer_size];

	int	a;
	int b;

	for(a=0 ; a!=4*mutation_Nmer_size ; ++a){
		new_ei_nucleotide_contributions[a] = ei_nucleotide_contributions[a];
	}

	for(a=0 ; a != mutation_Nmer_size ; ++a){
		//ei_nucleotide_contributions[i*mutation_Nmer_size+3] = 0;
		for(b=0 ; b!=3 ; ++b){
			//Update the contribution of the nucleotide
			new_ei_nucleotide_contributions[a*4+b] += alpha*gsl_vector_get(update_vect_p,(a*3 + b));

			//Compute the contribution of the constrained nucleotide
			new_ei_nucleotide_contributions[a*4+3] -= alpha*gsl_vector_get(update_vect_p,(a*3 + b));
		}
	}
	//Compute the values for the Jacobian and Hessian entries
	int base_4_address[mutation_Nmer_size];
	int max_address = 0;


	for (a=0;a!=mutation_Nmer_size;++a){
		base_4_address[a]=0;
		max_address += 3*adressing_vector[a];
	}
	max_address+=1;

	b=0;
	double current_Nmer_P_SHM;
	double current_Nmer_P_bg;
	double current_Nmer_unorm_score;
	double error_model_likelihood = 0;

	while(b!=max_address){
		//current_Nmer_P_SHM = Nmer_P_SHM[b];
		//current_Nmer_P_bg = Nmer_P_BG[b];
		current_Nmer_P_SHM = Nmer_N_SHM[b];
		current_Nmer_P_bg = Nmer_N_bg[b];
		current_Nmer_unorm_score = compute_Nmer_unorm_score(base_4_address,new_ei_nucleotide_contributions); //FIXME change this to chose the values of the eis and R to use


		//Update the base 10 and 4 addresses
		this->increment_base_10_and_4(b,base_4_address);

		error_model_likelihood+=current_Nmer_P_SHM*(log(new_R)+log(current_Nmer_unorm_score) - log(3)) - current_Nmer_P_bg*log(1+new_R*current_Nmer_unorm_score);
	}


	return error_model_likelihood;
}

void Hypermutation_global_errorrate::increment_base_10_and_4(int& base_10_counter , int* base_4_address){
	//Update the base 10 and 4 addresses
	++base_10_counter;//base 10
	bool bool_mod_4N = fmod(base_10_counter,4)==0;//FIXME define a function for this
	if(bool_mod_4N){
		int position = mutation_Nmer_size-1;
		while(bool_mod_4N){
			base_4_address[position]=0;
			--position;
			base_4_address[position]+=1;
			if(position==0)break;
			bool_mod_4N = fmod(base_10_counter,adressing_vector.at(position-1))==0;
			//Careful to out of range exception although should not happen
		}
	}
	else{
		base_4_address[mutation_Nmer_size-1]+=1;
	}
}

void Hypermutation_global_errorrate::initialize(const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map){
	//FIXME look for previous initialization to avoid memory leak

	//Initialize booleans for constructed sequences
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side))>0){
		v_gene=true;
	}
	else{v_gene=false;}
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side))>0){
		d_gene=true;
	}
	else{d_gene=false;}
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side))>0){
		j_gene=true;
	}
	else{j_gene=false;}
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Insertion_t,VJ_genes,Undefined_side))>0){
		vj_ins=true;
	}
	else{vj_ins=false;}
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Insertion_t,VD_genes,Undefined_side))>0){
		vd_ins=true;
	}
	else{vd_ins=false;}
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Insertion_t,DJ_genes,Undefined_side))>0){
		dj_ins=true;
	}
	else{dj_ins=false;}

	//Get the right pointers for the V gene
	if(v_gene){
		v_gene_event_p = dynamic_pointer_cast<Gene_choice> (events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side)));
		vgene_offset_p = &v_gene_event_p->alignment_offset_p;
		vgene_real_index_p = &v_gene_event_p->current_realization_index;

		//Initialize gene counters
		v_realizations = v_gene_event_p->get_realizations_map();
		//Get the number of realizations
		n_v_real = v_realizations.size();

		v_sequences = new Int_Str [n_v_real];
		for (const pair<const string,Event_realization> v_real: v_realizations){
			v_sequences[v_real.second.index] = v_real.second.value_str_int;
		}

		//Get deletion value pointer for V 3' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime)) != 0){
			shared_ptr<const Deletion> v_3_del_event_p = dynamic_pointer_cast<Deletion>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime)));
			v_3_del_value_p = &(v_3_del_event_p->deletion_value);
		}
		else{v_3_del_value_p = &no_del_buffer;}

	}
	else{
		if(learn_on == V_gene | learn_on == VJ_genes | learn_on == VD_genes | learn_on == VDJ_genes){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize V gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw runtime_error("Cannot learn on V gene without V choice in the model!");
		}
	}



	//Get the right pointers for the D gene
	if(d_gene){
		d_gene_event_p = dynamic_pointer_cast<Gene_choice>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side)));
		dgene_offset_p = &d_gene_event_p->alignment_offset_p;
		dgene_real_index_p = &d_gene_event_p->current_realization_index;

		//Get deletion value pointer for D 5' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime)) != 0){
			shared_ptr<const Deletion> d_5_del_event_p = dynamic_pointer_cast<Deletion>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime)));
			d_5_del_value_p = &(d_5_del_event_p->deletion_value);
		}
		else{d_5_del_value_p = &no_del_buffer;}

		//Get deletion value pointer for D 3' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime)) != 0){
			shared_ptr<const Deletion> d_3_del_event_p = dynamic_pointer_cast<Deletion>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime)));
			d_3_del_value_p = &(d_3_del_event_p->deletion_value);
		}
		else{d_3_del_value_p = &no_del_buffer;}

	}
	else{
		if(learn_on == D_gene | learn_on == DJ_genes | learn_on == VD_genes | learn_on == VDJ_genes){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize D gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw runtime_error("Cannot learn on D gene without D choice in the model!");
		}
	}


	//Get the right pointers for the J gene
	if(j_gene){
		j_gene_event_p = dynamic_pointer_cast<Gene_choice>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side)));
		jgene_offset_p = &j_gene_event_p->alignment_offset_p;
		jgene_real_index_p = &j_gene_event_p->current_realization_index;

		//Initialize gene counters
		j_realizations = j_gene_event_p->get_realizations_map();
		//Get the number of realizations
		n_j_real = j_realizations.size();

		j_sequences = new Int_Str [n_j_real];
		for (const pair<const string,Event_realization> j_real: j_realizations){
			j_sequences[j_real.second.index] = j_real.second.value_str_int;
		}

		//Get deletion value pointer for J 5' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime)) != 0){
			shared_ptr<const Deletion> j_5_del_event_p = dynamic_pointer_cast<Deletion>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime)));
			j_5_del_value_p = &(j_5_del_event_p->deletion_value);
		}
		else{j_5_del_value_p = &no_del_buffer;}
	}
	else{
		if(learn_on == J_gene | learn_on == DJ_genes | learn_on == VJ_genes | learn_on == VDJ_genes){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize J gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw runtime_error("Cannot learn on J gene without J choice in the model!");
		}
	}


	this->clean_all_counters();

}

void Hypermutation_global_errorrate::add_to_norm_counter(){
	if(seq_likelihood!=0){

		size_t array_size = pow(4,mutation_Nmer_size);
		for(size_t ii=0 ; ii != array_size ; ++ii){
			//cout<<debug_one_seq_Nmer_N_bg[ii]<<',';
			Nmer_N_SHM[ii] += one_seq_Nmer_N_SHM[ii]/seq_likelihood;
			one_seq_Nmer_N_SHM[ii]=0;

			Nmer_N_bg[ii] += one_seq_Nmer_N_bg[ii]/seq_likelihood;
			one_seq_Nmer_N_bg[ii]=0;
		}


		model_log_likelihood+=log10(seq_likelihood);
		number_seq+=1;


	}

	seq_mean_error_number = 0;
	seq_likelihood = 0;
	seq_probability = 0;
	debug_number_scenarios=0;


}

void Hypermutation_global_errorrate::clean_seq_counters(){
	//if(seq_likelihood!=0){
	size_t array_size = pow(4,mutation_Nmer_size);
	for(size_t ii=0 ; ii != array_size ; ++ii){
		//cout<<debug_one_seq_Nmer_N_bg[ii]<<',';
		one_seq_Nmer_N_SHM[ii]=0;

		one_seq_Nmer_N_bg[ii]=0;
	}



	seq_mean_error_number = 0;
	seq_likelihood = 0;
	seq_probability = 0;
	debug_number_scenarios=0;
}


void Hypermutation_global_errorrate::clean_all_counters(){

	size_t array_size = pow(4,mutation_Nmer_size);
	for(size_t ii=0 ; ii != array_size ; ++ii){
		Nmer_N_SHM[ii]=0;
		Nmer_N_bg[ii]=0;
	}

	this->clean_seq_counters();

}


void Hypermutation_global_errorrate::write2txt(ofstream& outfile){
	outfile<<"#Hypermutationglobalerrorrate;"<<this->mutation_Nmer_size<<";"<<this->learn_on<<";"<<this->apply_to<<endl;
	outfile<<mu<<endl;
	outfile<<ei_nucleotide_contributions[0];
	for(i=1 ; i!=mutation_Nmer_size*4 ; ++i){
		outfile<<";"<<ei_nucleotide_contributions[i];
	}
	outfile<<endl;
}

void Hypermutation_global_errorrate::update_Nmers_proba(int current_pos , int current_index,double current_score){
	//Iterate through possible nucleotides at this position
	for(int ii=0;ii!=4;++ii){
		int new_index = current_index + ii*adressing_vector[current_pos];
		double new_score = current_score*exp(ei_nucleotide_contributions[current_pos*4 +ii]);
		if(current_pos!=mutation_Nmer_size-1){
			update_Nmers_proba(current_pos+1,new_index,new_score);
		}
		else{
			this->Nmer_mutation_proba[new_index]= new_score*mu/(1+new_score*mu);
			//cout<<"Nmerproba: "<<new_index<<";"<<new_score*R/(1+new_score*R)<<endl;
		}
	}
}



/*void Hypermutation_global_errorrate::compute_P_SHM_and_BG(){
	//Initialize P_SHM and P_BG to 0


	if(learn_on_v){
		for(unordered_map<string,Event_realization>::const_iterator real_iter = v_realizations.begin() ; real_iter!=v_realizations.end() ; real_iter++){
			pair<size_t,double*> nucleotide_coverage = v_gene_nucleotide_coverage_p[(*real_iter).second.index];
			pair<size_t,double*> nucleotide_error = v_gene_per_nucleotide_error_p[(*real_iter).second.index];

			//Get the first Nmer on the gene
			Nmer_index = 0;
			current_Nmer = queue<size_t>();
			//and get min coverage for the first Nmer
			double min_coverage = INT16_MAX;
			for(j=0 ; j!=mutation_Nmer_size ; j++){
				tmp_int_nt = (*real_iter).second.value_str_int.at(j);
				current_Nmer.push(tmp_int_nt);
				Nmer_index+=adressing_vector[j]*tmp_int_nt;

				if(nucleotide_coverage.second[j]<min_coverage){
					min_coverage = nucleotide_coverage.second[j];
				}
			}

			Nmer_P_BG[Nmer_index] += min_coverage; //The coverage of the Nmer is only as high as the lowest covered nt
			Nmer_P_SHM[Nmer_index] += nucleotide_error.second[(mutation_Nmer_size-1)/2];

			for(i=(mutation_Nmer_size-1)/2 +1 ; i!=(*real_iter).second.value_str_int.size() - (mutation_Nmer_size-1)/2; ++i){

				//Remove the previous first nucleotide of the Nmer and it's contribution to the index
				Nmer_index-=current_Nmer.front()*adressing_vector[0];
				current_Nmer.pop();
				//Shift the index
				Nmer_index*=4;
				//Add the contribution of the new nucleotide

				tmp_int_nt = (*real_iter).second.value_str_int.at(i+(mutation_Nmer_size-1)/2);//Assume a symetric Nmer

				Nmer_index+=tmp_int_nt;
				current_Nmer.push(tmp_int_nt);

				double min_coverage = INT16_MAX;
				for(j=-(mutation_Nmer_size-1)/2 ; j!= (mutation_Nmer_size-1)/2 +1 ; j++){
					if(nucleotide_coverage.second[j+i]<min_coverage){
						min_coverage = nucleotide_coverage.second[j+i];
					}
				}

				Nmer_P_BG[Nmer_index] += min_coverage; //The coverage of the Nmer is only as high as the lowest covered nt
				Nmer_P_SHM[Nmer_index] += nucleotide_error.second[i];

			}
		}
	}

	if(learn_on_d){

	}
	//TODO remove code duplication this should be only one function taking nucleotide coverage and nucleotide error as arguments

	if(learn_on_j){
		for(unordered_map<string,Event_realization>::const_iterator real_iter = j_realizations.begin() ; real_iter!=j_realizations.end() ; real_iter++){
			pair<size_t,double*> nucleotide_coverage = j_gene_nucleotide_coverage_p[(*real_iter).second.index];
			pair<size_t,double*> nucleotide_error = j_gene_per_nucleotide_error_p[(*real_iter).second.index];

			//Get the first Nmer on the gene
			Nmer_index = 0;
			current_Nmer = queue<size_t>();
			//and get min coverage for the first Nmer
			double min_coverage = INT16_MAX;
			for(j=0 ; j!=mutation_Nmer_size ; j++){
				tmp_int_nt = (*real_iter).second.value_str_int.at(j);
				current_Nmer.push(tmp_int_nt);
				Nmer_index+=adressing_vector[j]*tmp_int_nt;

				if(nucleotide_coverage.second[j]<min_coverage){
					min_coverage = nucleotide_coverage.second[j];
				}
			}

			Nmer_P_BG[Nmer_index] += min_coverage; //The coverage of the Nmer is only as high as the lowest covered nt
			Nmer_P_SHM[Nmer_index] += nucleotide_error.second[(mutation_Nmer_size-1)/2];

			for(i=(mutation_Nmer_size-1)/2 +1 ; i!=(*real_iter).second.value_str_int.size() - (mutation_Nmer_size-1)/2; ++i){

				//Remove the previous first nucleotide of the Nmer and it's contribution to the index
				Nmer_index-=current_Nmer.front()*adressing_vector[0];
				current_Nmer.pop();
				//Shift the index
				Nmer_index*=4;
				//Add the contribution of the new nucleotide

				tmp_int_nt = (*real_iter).second.value_str_int.at(i+(mutation_Nmer_size-1)/2);//Assume a symetric Nmer

				Nmer_index+=tmp_int_nt;
				current_Nmer.push(tmp_int_nt);

				double min_coverage = INT16_MAX;
				for(j=-(mutation_Nmer_size-1)/2 ; j!= (mutation_Nmer_size-1)/2 +1 ; j++){
					if(nucleotide_coverage.second[j+i]<min_coverage){
						min_coverage = nucleotide_coverage.second[j+i];
					}
				}

				Nmer_P_BG[Nmer_index] += min_coverage; //The coverage of the Nmer is only as high as the lowest covered nt
				Nmer_P_SHM[Nmer_index] += nucleotide_error.second[i];

			}
		}
	}

	debug_stream<<endl<<"Pbg,Pshm,Nbg,Nshm"<<endl;
	for(int zzz=0 ; zzz!=pow(4,mutation_Nmer_size) ; ++zzz){
		debug_stream<<zzz<<","<<Nmer_P_BG[zzz]<<","<<Nmer_P_SHM[zzz]<<','<<debug_Nmer_N_bg[zzz]<<','<<debug_Nmer_N_SHM[zzz]<<endl;
	}
	debug_stream<<endl;
}*/



double Hypermutation_global_errorrate::compute_Nmer_unorm_score(int* base_4_address,double* ei_nucleotide_contributions_p ){
	double unorm_score = 1;
	for(int ii = 0 ; ii != mutation_Nmer_size ; ++ii){
		unorm_score*=exp(ei_nucleotide_contributions_p[4*ii+base_4_address[ii]]);
	}
	//cout<<"unorm score b4 address:"<<base_4_address[0]<<base_4_address[1]<<base_4_address[2]<<";"<<unorm_score<<endl;
	return unorm_score;
}

void Hypermutation_global_errorrate::introduce_uniform_transversion(char& nt , std::mt19937_64& generator , std::uniform_real_distribution<double>& distribution) const{
	double rand_trans = distribution(generator);

	if(nt == 'A'){
		if(rand_trans<= 1.0/3.0){
			nt = 'C';
		}
		else if (rand_trans >= 2.0/3.0){
			nt = 'G';
		}
		else{
			nt = 'T';
		}
	}
	else if(nt == 'C'){
		if(rand_trans<= 1.0/3.0){
			nt = 'A';
		}
		else if (rand_trans >= 2.0/3.0){
			nt = 'G';
		}
		else{
			nt = 'T';
		}
	}
	else if(nt == 'G'){
		if(rand_trans<= 1.0/3.0){
			nt = 'A';
		}
		else if (rand_trans >= 2.0/3.0){
			nt = 'C';
		}
		else{
			nt = 'T';
		}

	}
	else if (nt == 'T'){
		if(rand_trans<= 1.0/3.0){
			nt = 'A';
		}
		else if (rand_trans >= 2.0/3.0){
			nt = 'C';
		}
		else{
			nt = 'G';
		}
	}
	else{
		throw runtime_error("unknown nucleotide in Hypermutationglobalerrorrate::generate_errors()");
	}
}
