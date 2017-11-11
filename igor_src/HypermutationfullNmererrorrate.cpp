/*
 * HypermutationfullNmererrorrate.cpp
 *
 *  Created on: Aug 30, 2017
 *      Author: quentin
 */

#include "HypermutationfullNmererrorrate.h"

using namespace std;

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value,size_t n_observed_thresh/*=0*/): Error_rate() , n_observed_Nmer_threshold(n_observed_thresh) , mutation_Nmer_size(nmer_width) , learn_on(learn) , apply_to(apply) , n_v_real(0) , n_j_real(0) , n_d_real(0) ,
		v_sequences(NULL),j_sequences(NULL),
		dj_ins(true) , vd_ins(true) , vj_ins(true) , v_gene(true) , d_gene(true) , j_gene(true) ,
		vgene_offset_p(NULL) , dgene_offset_p(NULL) , jgene_offset_p(NULL) ,
		vgene_real_index_p(NULL) , dgene_real_index_p(NULL) , jgene_real_index_p(NULL),
		v_3_del_value_p(NULL) , d_5_del_value_p(NULL) , d_3_del_value_p(NULL) , j_5_del_value_p(NULL),
		i(-1) , j(-1) , v_3_del_value_corr(INT16_MAX) , d_5_del_value_corr(INT16_MAX) , d_3_del_value_corr(INT16_MAX) , j_5_del_value_corr(INT16_MAX) , tmp_cov_p(NULL) , tmp_err_p(NULL) , tmp_corr_len(-1) , tmp_len_util(-1) , scenario_new_proba(-1) ,
		largest_nuc_adress(-1), tmp_int_nt(-1) , Nmer_index(-1),
		output_Nmer_stat_stream(new ofstream){


	if(fmod(nmer_width,2)==0){
		throw runtime_error("Cannot instanciate hypermutation full Nmer error rate with an even size Nmer(need to be symmetric) in Hypermutation_full_Nmer_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value)");
	}

	size_t array_size = pow(4,mutation_Nmer_size);

	if( (starting_flat_value<0) or (starting_flat_value>1) ){
		throw invalid_argument("The starting flat value for the hypermutation probability must lie between 0 and 1, passed value is " + to_string(starting_flat_value) + " in Hypermutation_full_Nmer_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value)");
	}

	// Instantiate and initialize arrays
	Nmer_mutation_proba = new double [array_size];
	one_seq_Nmer_N_SHM = new double [array_size];
	one_seq_Nmer_N_bg = new double [array_size];
	Nmer_N_SHM = new double [array_size];
	Nmer_N_bg = new double [array_size];
	for(int ii=0 ; ii!=array_size ; ++ii){
		Nmer_mutation_proba[ii] = starting_flat_value;
		one_seq_Nmer_N_SHM[ii] = 0;
		one_seq_Nmer_N_bg[ii] = 0;
		Nmer_N_SHM[ii] = 0;
		Nmer_N_bg[ii] = 0;
	}

	//Now the the probability array is initialized, build the upper bound matrix
	build_upper_bound_matrix(1,1);


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

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , vector<double> init_Nmer_mutations_probas,size_t n_observed_thresh/*=0*/): Hypermutation_full_Nmer_errorrate(nmer_width , learn , apply , 0,n_observed_thresh){
	if(init_Nmer_mutations_probas.size()==pow(4,mutation_Nmer_size)){
		for(i=0 ; i != init_Nmer_mutations_probas.size() ; ++i){
			if((init_Nmer_mutations_probas[i]>=0) and (init_Nmer_mutations_probas[i]<=1)){
				this->Nmer_mutation_proba[i] = init_Nmer_mutations_probas[i];
			}
			else{
				throw invalid_argument("The starting values for the hypermutation probabilities must lie between 0 and 1, passed value is " + to_string(init_Nmer_mutations_probas[i]) + "for Nmer index " + to_string(i) + " in Hypermutation_full_Nmer_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , vector<double> init_Nmer_mutations_probas)");
			}
		}
	}
	else{
		throw runtime_error("Size of Nmer mutation probabilities vector does not match the expected size in Hypermutation_full_Nmer_errorrate(size_t,Gene_class,Gene_class,double,std::vector<double>)");
	}
}

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , double starting_flat_value,string filename,size_t n_observed_thresh/*=0*/):Hypermutation_full_Nmer_errorrate(nmer_width , learn , apply , starting_flat_value,n_observed_thresh) {
	this->set_output_Nmer_stream(filename);
}

Hypermutation_full_Nmer_errorrate::Hypermutation_full_Nmer_errorrate(size_t nmer_width , Gene_class learn , Gene_class apply , vector<double> init_Nmer_mutations_probas,string filename,size_t n_observed_thresh/*=0*/):Hypermutation_full_Nmer_errorrate(nmer_width , learn , apply ,init_Nmer_mutations_probas,n_observed_thresh) {
	this->set_output_Nmer_stream(filename);
}

Hypermutation_full_Nmer_errorrate::~Hypermutation_full_Nmer_errorrate() {
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

void Hypermutation_full_Nmer_errorrate::set_output_Nmer_stream(string filename){
	cout<<"Full Nmer hypermutation model output set to: "<<filename<<endl;
	output_Nmer_stat_stream->open(filename);
	(*output_Nmer_stat_stream)<<"Nmer_index,N_bg,N_mut"<<endl;
	output_Nmer_stat = true;
}

shared_ptr<Error_rate> Hypermutation_full_Nmer_errorrate::copy()const{

	shared_ptr<Hypermutation_full_Nmer_errorrate> copy_err_r = shared_ptr<Hypermutation_full_Nmer_errorrate>( new Hypermutation_full_Nmer_errorrate(this->mutation_Nmer_size , this->learn_on , this->apply_to , 0 , this->n_observed_Nmer_threshold) );
	copy_err_r->updated = this->updated;
	copy_err_r->output_Nmer_stat = this->output_Nmer_stat;
	copy_err_r->output_Nmer_stat_stream = this->output_Nmer_stat_stream;
	//copy_err_r->R = this->R;
	for(int ii = 0 ; ii != pow(4,mutation_Nmer_size) ; ++ii){
		copy_err_r->Nmer_mutation_proba[ii] = this->Nmer_mutation_proba[ii];
	}

	return copy_err_r;

}

Hypermutation_full_Nmer_errorrate& Hypermutation_full_Nmer_errorrate::operator +=(Hypermutation_full_Nmer_errorrate err_r){

	// CHeck whether all mutations probas are the same
	bool identical_mut_probas = true;
	for(int ii = 0 ; ii != pow(4,mutation_Nmer_size) ; ++ii){
		if(err_r.Nmer_mutation_proba[ii] != this->Nmer_mutation_proba[ii]){
			identical_mut_probas = false;
		}
	}

	//FIXME sequential ifs throwing more meaningful exception
	if( (this->learn_on == err_r.learn_on)
		& (this->apply_to == err_r.apply_to)
		& (this->mutation_Nmer_size == err_r.mutation_Nmer_size)
		& identical_mut_probas){
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
		throw runtime_error("Hypermutation models cannot be added in Hypermutation_full_Nmer_errorrate::operator +=()");
	}

}

Error_rate* Hypermutation_full_Nmer_errorrate::add_checked(Error_rate* err_r){
	return &(this->operator +=( *(dynamic_cast<Hypermutation_full_Nmer_errorrate*>(err_r) ) ));
}

const double& Hypermutation_full_Nmer_errorrate::get_err_rate_upper_bound(size_t n_errors , size_t n_error_free) {


	if( n_errors>this->max_err || n_error_free>this->max_noerr){
		/*double max_mut_proba = 0;
		double min_mut_proba = 1;
		size_t array_size = pow(4,mutation_Nmer_size);
		for(size_t ii=0 ; ii != array_size ; ++ii){
			if(max_mut_proba<this->Nmer_mutation_proba[ii]) max_mut_proba=this->Nmer_mutation_proba[ii];
			if(min_mut_proba>this->Nmer_mutation_proba[ii]) min_mut_proba=this->Nmer_mutation_proba[ii];
		}*/

		size_t array_size = pow(4,mutation_Nmer_size);
		vector<double> probas_vector;
		for(size_t ii=0 ; ii != array_size ; ++ii){
			probas_vector.push_back(this->Nmer_mutation_proba[ii]);
		}
		sort(probas_vector.begin(),probas_vector.end());
		//By definition the number of mutation probabilities is even (power of 4)
		double median_mut_proba = (probas_vector[probas_vector.size()/2 -1] + probas_vector[probas_vector.size()/2])/2.0;

		//Need to increase the matrix size (anyway the matrix is at very most read_len^2
		Matrix<double> new_bound_mat (max(this->max_err,n_errors + 10) , max(this->max_noerr , n_error_free+10));
		for(size_t i=0 ; i!=new_bound_mat.get_n_rows() ; ++i){
			for(size_t j=0 ; j!=new_bound_mat.get_n_cols() ; ++j){
				if(i<this->max_err and j<this->max_noerr){
					new_bound_mat(i,j) = this->upper_bound_proba_mat(i,j);
				}
				else{
					new_bound_mat(i,j) = pow(median_mut_proba/3.0,i)*pow((1-median_mut_proba),j);
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

void Hypermutation_full_Nmer_errorrate::build_upper_bound_matrix(size_t m,size_t n){
	Matrix<double> new_bound_mat (m,n);

	//Get min and max mutation proba
	/*double max_mut_proba = 0;
	double min_mut_proba = 1;
	size_t array_size = pow(4,mutation_Nmer_size);
	for(size_t ii=0 ; ii != array_size ; ++ii){
		if(max_mut_proba<this->Nmer_mutation_proba[ii]) max_mut_proba=this->Nmer_mutation_proba[ii];
		if(min_mut_proba>this->Nmer_mutation_proba[ii]) min_mut_proba=this->Nmer_mutation_proba[ii];
	}*/

	//Get the median hypermutation probability
	size_t array_size = pow(4,mutation_Nmer_size);
	vector<double> probas_vector;
	for(size_t ii=0 ; ii != array_size ; ++ii){
		probas_vector.push_back(this->Nmer_mutation_proba[ii]);
	}
	sort(probas_vector.begin(),probas_vector.end());
	//By definition the number of mutation probabilities is even (power of 4)
	double median_mut_proba = (probas_vector[probas_vector.size()/2 -1] + probas_vector[probas_vector.size()/2])/2.0;

	//Fill the matrix
	for(size_t i=0 ; i!=new_bound_mat.get_n_rows() ; ++i){
		for(size_t j=0 ; j!=new_bound_mat.get_n_cols() ; ++j){
			if(i<this->max_err and j<this->max_noerr){
				new_bound_mat(i,j) = this->upper_bound_proba_mat(i,j);
			}
			else{
				//new_bound_mat(i,j) = pow(max_mut_proba/3.0,i)*pow((1-min_mut_proba),j);
				new_bound_mat(i,j) = pow(median_mut_proba/3.0,i)*pow((1-median_mut_proba),j);
			}
		}
	}
	this->upper_bound_proba_mat = new_bound_mat;
	this->max_err = new_bound_mat.get_n_rows()-1;
	this->max_noerr = new_bound_mat.get_n_cols()-1;
}

double Hypermutation_full_Nmer_errorrate::compare_sequences_error_prob (double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){
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
			tmp_len_util = seq_offsets.at(J_gene_seq,Five_prime)-(mutation_Nmer_size-1)/2; //Start using the information of the (N-1)/2 inserted (or D) nucleotides before the J

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

	}


	this->seq_likelihood += scenario_new_proba;
	this->seq_probability+=scenario_probability;
	this->seq_mean_error_number +=  (v_mismatch_list.size() + d_mismatch_list.size() + j_mismatch_list.size())*scenario_new_proba;
	++debug_number_scenarios;

	return scenario_new_proba;

}

queue<int> Hypermutation_full_Nmer_errorrate::generate_errors(string& generated_seq , default_random_engine& generator) const{
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

unsigned Hypermutation_full_Nmer_errorrate::generate_random_mutation_probas(double mean, double std){
	//Create seed for random generator
	//create a seed from timer
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point time = myclock::now();
	myclock::duration dur = myclock::time_point::max() - time;

	unsigned time_seed = dur.count();
	//Instantiate random number generator
	default_random_engine generator =  default_random_engine(time_seed);
	normal_distribution<double> distribution(mean,std);

	size_t array_size = pow(4,mutation_Nmer_size);
	for(i = 0 ; i != array_size ; ++i){
		Nmer_mutation_proba[i] =  distribution(generator);
	}

	return time_seed;
}


void Hypermutation_full_Nmer_errorrate::update(){
	size_t array_size = pow(4,mutation_Nmer_size);
	//If an output stream is provided outputs Nmer statistics in a file (mostly for debugging)
	if(output_Nmer_stat){
		for(size_t zzz=0 ; zzz!=array_size ; ++zzz){
			(*output_Nmer_stat_stream)<<zzz<<","<<Nmer_N_bg[zzz]<<","<<Nmer_N_SHM[zzz]<<endl;
		}
	}


	// Update the error rate by maximizing the likelihood of the error model
	// This simply boils down to equating the model mutation probabilities to the posterior mutation frequencies

	//double average_mutability = 0;
	//size_t n_observed_nmers = 0;
	vector<double> trusted_probas_vector;
	for(size_t ii = 0 ; ii!=array_size ; ++ii){
		//Only update the value if the Nmer has been observed
		//Note that if an Nmer is not observed much the mutation probability might artificially go to 0 because of undersampling, remain cautious when interpreting such values
		if(Nmer_N_bg[ii]>=n_observed_Nmer_threshold){
			Nmer_mutation_proba[ii] = Nmer_N_SHM[ii]/Nmer_N_bg[ii];
			trusted_probas_vector.push_back(this->Nmer_mutation_proba[ii]); //Compute the median mutability value over trustworthy Nmers
			//average_mutability+=Nmer_mutation_proba[ii];
			//++n_observed_nmers;
		}
	}
	sort(trusted_probas_vector.begin(),trusted_probas_vector.end());
	//Get the median
	double median_mut_proba;
	if( (trusted_probas_vector.size()%2) == 0){
		median_mut_proba = (trusted_probas_vector[trusted_probas_vector.size()/2 -1] + trusted_probas_vector[trusted_probas_vector.size()/2])/2.0;
	}
	else{
		median_mut_proba = trusted_probas_vector[trusted_probas_vector.size()/2];
	}



	//Now replace the value by the average mutability for unobserved Nmers
	//average_mutability/=n_observed_nmers;
	for(size_t ii = 0 ; ii!=array_size ; ++ii){
		if(Nmer_N_bg[ii]<n_observed_Nmer_threshold){
			Nmer_mutation_proba[ii] = median_mut_proba;
		}
	}

	//Clean counters
	this->clean_all_counters();

}


void Hypermutation_full_Nmer_errorrate::initialize(const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map){
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

void Hypermutation_full_Nmer_errorrate::add_to_norm_counter(){
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

void Hypermutation_full_Nmer_errorrate::clean_seq_counters(){
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


void Hypermutation_full_Nmer_errorrate::clean_all_counters(){

	size_t array_size = pow(4,mutation_Nmer_size);
	for(size_t ii=0 ; ii != array_size ; ++ii){
		Nmer_N_SHM[ii]=0;
		Nmer_N_bg[ii]=0;
	}

	this->clean_seq_counters();

}


void Hypermutation_full_Nmer_errorrate::write2txt(ofstream& outfile){
	outfile<<"#HypermutationfullNmererrorrate;"<<this->mutation_Nmer_size<<";"<<this->learn_on<<";"<<this->apply_to<<endl;
	outfile<<Nmer_mutation_proba[0];
	for(i=1 ; i!=pow(4,mutation_Nmer_size) ; ++i){
		outfile<<";"<<Nmer_mutation_proba[i];
	}
	outfile<<endl;
}

void Hypermutation_full_Nmer_errorrate::introduce_uniform_transversion(char& nt , std::default_random_engine& generator , std::uniform_real_distribution<double>& distribution) const{
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

