/*
 * Singleerrorrate.cpp
 *
 *  Created on: Jan 23, 2015
 *      Author: quentin
 */

#include "Singleerrorrate.h"

using namespace std;

Single_error_rate::Single_error_rate(): model_rate(0) , normalized_counter(0) , number_seq(0) , seq_weighted_er(0) {}

Single_error_rate::Single_error_rate(double error_rate): model_rate(error_rate) , normalized_counter(0) , number_seq(0) , seq_weighted_er(0)  {}

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

double Single_error_rate::get_err_rate_upper_bound() const{
	return this->model_rate/3;
	//TODO remove factor 3
}


double Single_error_rate::compare_sequences_error_prob (double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){
	//TODO extract sequence comparision from here, implement it in Errorrate class??
	number_errors=0;
	//cout<<constructed_sequences.at(V_gene_seq);
	genomic_nucl=0;

	string& v_gene_seq = (*constructed_sequences[V_gene_seq]);
	string& d_gene_seq = (*constructed_sequences[D_gene_seq]);
	string& j_gene_seq = (*constructed_sequences[J_gene_seq]);

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
		this->seq_weighted_er += temp;
		this->seq_likelihood += scenario_new_proba;
		this->seq_probability+=scenario_probability;
		++debug_number_scenarios;
		return scenario_new_proba;
	}
	else{
		return 0;
	}

}

queue<int> Single_error_rate::generate_errors(string& generated_seq , default_random_engine& generator) const{
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
