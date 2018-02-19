/*
 * Singleerrorrate.cpp
 *
 *  Created on: Jan 23, 2015
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

#include "Singleerrorrate.h"

using namespace std;

Single_error_rate::Single_error_rate(): Single_error_rate(0.0) {}

Single_error_rate::Single_error_rate(double error_rate): Error_rate() , model_rate(error_rate) , normalized_counter(0) , seq_weighted_er(0)  {
	build_upper_bound_matrix(1,1);
}

Single_error_rate::~Single_error_rate() {
	// TODO Auto-generated destructor stub
}

Single_error_rate Single_error_rate::operator +(Single_error_rate err_r){
	Single_error_rate temp =  *this;
	return temp+=err_r;

}

Single_error_rate& Single_error_rate::operator +=(Single_error_rate err_r){
	this->normalized_counter+= err_r.normalized_counter;
	this->number_seq+=err_r.number_seq;
	this->model_log_likelihood+=err_r.model_log_likelihood;
	return *this;
}

shared_ptr<Error_rate> Single_error_rate::copy()const{
	shared_ptr<Single_error_rate> copy_err_r = shared_ptr<Single_error_rate>(new Single_error_rate(this->model_rate));
	copy_err_r->updated = this->updated;
	return copy_err_r;
}

Error_rate* Single_error_rate::add_checked(Error_rate* err_r){
	return &( this->operator +=( *( dynamic_cast< Single_error_rate*> (err_r) ) ) );
}

const double& Single_error_rate::get_err_rate_upper_bound(size_t n_errors , size_t n_error_free) {
	if( n_errors>this->max_err || n_error_free>this->max_noerr){
		//Need to increase the matrix size (anyway the matrix is at very most read_len^2
		this->build_upper_bound_matrix(max(this->max_err,n_errors + 10) , max(this->max_noerr , n_error_free+10));
	}

	return this->upper_bound_proba_mat(n_errors,n_error_free);

	//return this->model_rate/3;
	//TODO remove factor 3
}

void Single_error_rate::build_upper_bound_matrix(size_t m , size_t n){
	Matrix<double> new_bound_mat (m,n);
	for(size_t i=0 ; i!=new_bound_mat.get_n_rows() ; ++i){
		for(size_t j=0 ; j!=new_bound_mat.get_n_cols() ; ++j){
			if(i<this->max_err and j<this->max_noerr){
				new_bound_mat(i,j) = this->upper_bound_proba_mat(i,j);
			}
			else{
				new_bound_mat(i,j) = pow(this->model_rate/3,i)*pow(1-this->model_rate,j);
			}
			//new_bound_mat(i,j) = pow(this->model_rate/3,i)*pow(1-this->model_rate,j);
		}
	}
	this->upper_bound_proba_mat = new_bound_mat;
	this->max_err = new_bound_mat.get_n_rows()-1;
	this->max_noerr = new_bound_mat.get_n_cols()-1;
}


double Single_error_rate::compare_sequences_error_prob (double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){
	//TODO extract sequence comparision from here, implement it in Errorrate class??
	number_errors=0;
	//cout<<constructed_sequences.at(V_gene_seq);
	genomic_nucl=0;

	Int_Str& v_gene_seq = (*constructed_sequences[V_gene_seq]);
	Int_Str& d_gene_seq = (*constructed_sequences[D_gene_seq]);
	Int_Str& j_gene_seq = (*constructed_sequences[J_gene_seq]);

	vector<int>& v_mismatch_list = *mismatches_lists[V_gene_seq];
	if(mismatches_lists.exist(D_gene_seq)){
		vector<int>& d_mismatch_list = *mismatches_lists[D_gene_seq];
		number_errors+=d_mismatch_list.size();
		genomic_nucl+=d_gene_seq.size();
	}

	vector<int>& j_mismatch_list = *mismatches_lists[J_gene_seq];

	genomic_nucl+=v_gene_seq.size();
	//genomic_nucl+=d_gene_seq.size();
	genomic_nucl+=j_gene_seq.size();

	number_errors+=v_mismatch_list.size();
	//number_errors+=d_mismatch_list.size();
	number_errors+=j_mismatch_list.size();


	 // Here a long double is required in case a lot of errors occur and/or the model rate is low, the probability will be truncated to 0 if it gets below ± 2.225,073,858,507,201,4 · 10-308 with double precision



	scenario_new_proba = scenario_probability*pow(model_rate/3,number_errors)*pow(1-model_rate,genomic_nucl-number_errors);
	if(scenario_new_proba >= seq_max_prob_scenario*proba_threshold_factor) {
		//if genomic nucl != 0 ?
		this->seq_mean_error_number +=  number_errors*scenario_new_proba;
		temp2 = (double(number_errors)/double(genomic_nucl));
		temp = scenario_new_proba*temp2;
		if(viterbi_run){
			this->seq_weighted_er = temp;
			this->seq_likelihood = scenario_new_proba;
			this->seq_probability = scenario_probability;
		}
		else{
			this->seq_weighted_er += temp;
			this->seq_likelihood += scenario_new_proba;
			this->seq_probability += scenario_probability;
		}

		++debug_number_scenarios;
		return scenario_new_proba;
	}
	else{
		return 0;
	}

}

queue<int> Single_error_rate::generate_errors(string& generated_seq , mt19937_64& generator) const{
	uniform_real_distribution<double> distribution(0.0,1.0);
	double rand_err ;// distribution(generator);
	double rand_trans;
	size_t index = 0;
	queue<int> errors_indices;
	for(string::iterator iter = generated_seq.begin() ; iter != generated_seq.end() ; ++iter){
		rand_err = distribution(generator);
		if(rand_err<=this->model_rate){
			//Introduce an error
			rand_trans = distribution(generator);
			errors_indices.push(index);

			if((*iter) == 'A'){
				if(rand_trans<= 1.0/3.0){
					(*iter) = 'C';
				}
				else if (rand_trans >= 2.0/3.0){
					(*iter) = 'G';
				}
				else{
					(*iter) = 'T';
				}
			}
			else if((*iter) == 'C'){
				if(rand_trans<= 1.0/3.0){
					(*iter) = 'A';
				}
				else if (rand_trans >= 2.0/3.0){
					(*iter) = 'G';
				}
				else{
					(*iter) = 'T';
				}
			}
			else if((*iter) == 'G'){
				if(rand_trans<= 1.0/3.0){
					(*iter) = 'A';
				}
				else if (rand_trans >= 2.0/3.0){
					(*iter) = 'C';
				}
				else{
					(*iter) = 'T';
				}

			}
			else if ((*iter == 'T')){
				if(rand_trans<= 1.0/3.0){
					(*iter) = 'A';
				}
				else if (rand_trans >= 2.0/3.0){
					(*iter) = 'C';
				}
				else{
					(*iter) = 'G';
				}

			}
			else{
				throw runtime_error("unknown nucleotide in Single_error_rate::generate_errors()");
			}

		}
		else{
			//Do nothing
		}
		++index;
	}
	return errors_indices;
}

int Single_error_rate::subseq_compare_err_num(const string& original_sequence , const string& constructed_sequence){
	int number_errors=0;

	for (size_t nucl_ind = 0 ; nucl_ind != original_sequence.size() ; ++nucl_ind){
			//Only take into account genomic nucl
			//might have to refine this if  the type of the inserted nucleotide is inferred (use constructed sequences map)
			if(original_sequence.at(nucl_ind) != constructed_sequence.at(nucl_ind)){
				number_errors+=1;
			}
		}
	return number_errors;
}

void Single_error_rate::update(){
	if(this->is_updated()){
		model_rate = normalized_counter / number_seq;
		normalized_counter = 0;
		number_seq = 0;
	}
}
/*
 * This method ensure correct normalization of the error rate
 * The error rate inferred for each scenario is weighted by the probability of the scenario.
 * The sum of these is then normalized by the likelihood of the sequence
 */
void Single_error_rate::add_to_norm_counter(){

	if(seq_likelihood != 0){ //TODO check that the first version was not more correct
			normalized_counter += seq_weighted_er/seq_likelihood;
			model_log_likelihood+=log10(seq_likelihood);
			number_seq += 1;
	}

	/*if(seq_weighted_er != 0){ //Why seq weighted error instead of seq likelihood?
		normalized_counter += seq_weighted_er/seq_likelihood;
		model_log_likelihood+=log10(seq_likelihood);
	}*/
	seq_weighted_er = 0;
	seq_likelihood = 0;
	seq_probability = 0;
	debug_number_scenarios=0;
	seq_mean_error_number=0;

}
/*
 * This method cleans the sequence specific counters
 * This method is used in case the sequence has a mean number of errors greater than the threshold
 * In this case the sequence will not contribute to the error rate
 */
void Single_error_rate::clean_seq_counters(){
	seq_weighted_er = 0;
	seq_likelihood = 0;
	seq_probability = 0;
	debug_number_scenarios=0;
	seq_mean_error_number=0;
}

void Single_error_rate::write2txt(ofstream& outfile){
	outfile<<"#SingleErrorRate"<<endl;
	outfile<<model_rate<<endl;
}
