/*
 * Coverageerrcounter.cpp
 *
 *  Created on: Oct 24, 2016
 *      Author: quentin
 */

#include "Coverageerrcounter.h"

using namespace std;


Coverage_err_counter::Coverage_err_counter(Gene_class count_on): Coverage_err_counter("/tmp/",count_on,1,false,false){

}

Coverage_err_counter::Coverage_err_counter(Gene_class count_on , bool dump_all_seq , bool last_iter_only): Coverage_err_counter("/tmp/",count_on,1,dump_all_seq,last_iter_only){

}


Coverage_err_counter::Coverage_err_counter(string path , Gene_class count_on , bool dump_all_seq): Coverage_err_counter(path,count_on,1,dump_all_seq,false){

}

Coverage_err_counter::Coverage_err_counter(string path , Gene_class count_on , size_t Npoint_count , bool dump_all_seq , bool last_iter_only): Counter(path , last_iter_only) ,
		count_on(count_on) , dump_individual_seqs(dump_all_seq), record_Npoint_occurence(Npoint_count),
		positions(NULL),
		n_v_real(0),n_d_real(0),n_j_real(0),
		v_gene_nucleotide_coverage_p(NULL) , v_gene_per_nucleotide_error_p(NULL),d_gene_nucleotide_coverage_p(NULL) , d_gene_per_nucleotide_error_p(NULL),j_gene_nucleotide_coverage_p(NULL) , j_gene_per_nucleotide_error_p(NULL),
		v_gene_nucleotide_coverage_seq_p(NULL) , v_gene_per_nucleotide_error_seq_p(NULL) , d_gene_nucleotide_coverage_seq_p(NULL) , d_gene_per_nucleotide_error_seq_p(NULL) , j_gene_nucleotide_coverage_seq_p(NULL) , j_gene_per_nucleotide_error_seq_p(NULL) ,
		vgene_offset_p(NULL) , dgene_offset_p(NULL) , jgene_offset_p(NULL) ,
		vgene_real_index_p(NULL) , dgene_real_index_p(NULL) , jgene_real_index_p(NULL),
		v_3_del_value_p(NULL) , d_5_del_value_p(NULL) , d_3_del_value_p(NULL) , j_5_del_value_p(NULL){

	if(count_on == V_gene | count_on == VJ_genes | count_on == VD_genes | count_on == VDJ_genes){
			count_on_v = true;
		}
		else count_on_v = false;

		if(count_on == D_gene | count_on == DJ_genes | count_on == VD_genes | count_on == VDJ_genes){
			count_on_d = true;
		}
		else count_on_d = false;

		if(count_on == J_gene | count_on == VJ_genes | count_on == DJ_genes | count_on == VDJ_genes){
			count_on_j = true;
		}
		else count_on_j = false;

}

Coverage_err_counter::~Coverage_err_counter() {
	//FIXME TONS OF STUFF TO DELETE
	if(count_on_v){
		this->deallocate_coverage_and_errors_arrays(n_v_real,v_realizations,v_gene_nucleotide_coverage_p,v_gene_per_nucleotide_error_p,v_gene_nucleotide_coverage_seq_p,v_gene_per_nucleotide_error_seq_p);
	}
	if(count_on_d){
		this->deallocate_coverage_and_errors_arrays(n_d_real,d_realizations,d_gene_nucleotide_coverage_p,d_gene_per_nucleotide_error_p,d_gene_nucleotide_coverage_seq_p,d_gene_per_nucleotide_error_seq_p);
	}
	if(count_on_j){
		this->deallocate_coverage_and_errors_arrays(n_j_real,j_realizations,j_gene_nucleotide_coverage_p,j_gene_per_nucleotide_error_p,j_gene_nucleotide_coverage_seq_p,j_gene_per_nucleotide_error_seq_p);
	}

}

void Coverage_err_counter::initialize_counter(const Model_Parms& parms , const Model_marginals& marginals){
	if(not fstreams_created){

		if(count_on_v){
			output_cov_err_v_file_ptr = shared_ptr<ofstream>(new ofstream);
			this->output_cov_err_v_file_ptr->open(path_to_file + "V_genes_cov_and_err.csv");

			//Create the header
			if(dump_individual_seqs){
				(*this->output_cov_err_v_file_ptr.get())<<"iteration_n;seq_index;gene_index;coverage;errors"<<endl;
			}
			else{
				(*this->output_cov_err_v_file_ptr.get())<<"iteration_n;gene_index;coverage;errors"<<endl;
			}
		}

		if(count_on_d){
			output_cov_err_d_file_ptr = shared_ptr<ofstream>(new ofstream);
			this->output_cov_err_d_file_ptr->open(path_to_file + "D_genes_cov_and_err.csv");

			//Create the header
			if(dump_individual_seqs){
				(*this->output_cov_err_d_file_ptr.get())<<"iteration_n;seq_index;gene_index;coverage;errors"<<endl;
			}
			else{
				(*this->output_cov_err_d_file_ptr.get())<<"iteration_n;gene_index;coverage;errors"<<endl;
			}
		}

		if(count_on_j){
			output_cov_err_j_file_ptr = shared_ptr<ofstream>(new ofstream);
			this->output_cov_err_j_file_ptr->open(path_to_file + "J_genes_cov_and_err.csv");

			//Create the header
			if(dump_individual_seqs){
				(*this->output_cov_err_j_file_ptr.get())<<"iteration_n;seq_index;gene_index;coverage;errors"<<endl;
			}
			else{
				(*this->output_cov_err_j_file_ptr.get())<<"iteration_n;gene_index;coverage;errors"<<endl;
			}
		}

		fstreams_created = true;

	}
	positions = new size_t [record_Npoint_occurence];
	auto events_map = parms.get_events_map();

	if(count_on_v){
		//Initialize V pointers
		try{
			v_gene_event_p = dynamic_pointer_cast<Gene_choice> (events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side)));
			vgene_offset_p = &v_gene_event_p->alignment_offset_p;
			vgene_real_index_p = &v_gene_event_p->current_realization_index;

			//Initialize gene counters
			v_realizations = v_gene_event_p->get_realizations_map();
			//Get the number of realizations
			n_v_real = v_realizations.size();


			this->allocate_coverage_and_errors_arrays(n_v_real,v_realizations,v_gene_nucleotide_coverage_p,v_gene_per_nucleotide_error_p,v_gene_nucleotide_coverage_seq_p,v_gene_per_nucleotide_error_seq_p);

		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize V gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}

		//Get deletion value pointer for V 3' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime)) != 0){
			shared_ptr<const Deletion> v_3_del_event_p = dynamic_pointer_cast<Deletion>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime)));
			v_3_del_value_p = &(v_3_del_event_p->deletion_value);
		}
		else{v_3_del_value_p = &no_del_buffer;}

	}


	if(count_on_d){
		//Initialize D pointers
		try{
			d_gene_event_p = dynamic_pointer_cast<Gene_choice> (events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side)));
			dgene_offset_p = &d_gene_event_p->alignment_offset_p;
			dgene_real_index_p = &d_gene_event_p->current_realization_index;

			//Initialize gene counters
			d_realizations = d_gene_event_p->get_realizations_map();
			//Get the number of realizations
			n_d_real = d_realizations.size();


			this->allocate_coverage_and_errors_arrays(n_d_real,d_realizations,d_gene_nucleotide_coverage_p,d_gene_per_nucleotide_error_p,d_gene_nucleotide_coverage_seq_p,d_gene_per_nucleotide_error_seq_p);

		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize D gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}


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

	if(count_on_j){
		//Initialize J pointers
		try{
			j_gene_event_p = dynamic_pointer_cast<Gene_choice> (events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side)));
			jgene_offset_p = &j_gene_event_p->alignment_offset_p;
			jgene_real_index_p = &j_gene_event_p->current_realization_index;

			//Initialize gene counters
			j_realizations = j_gene_event_p->get_realizations_map();
			//Get the number of realizations
			n_j_real = j_realizations.size();


			this->allocate_coverage_and_errors_arrays(n_j_real,j_realizations,j_gene_nucleotide_coverage_p,j_gene_per_nucleotide_error_p,j_gene_nucleotide_coverage_seq_p,j_gene_per_nucleotide_error_seq_p);

		}
		catch(exception& except){
			cout<<"Exception caught during initialization of Hypermutation global error rate"<<endl;
			cout<<"Exception caught trying to initialize J gene pointers"<<endl;
			cout<<endl<<"throwing exception now..."<<endl;
			throw except;
		}

		//Get deletion value pointer for J 5' deletions if it exists
		if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime)) != 0){
			shared_ptr<const Deletion> j_5_del_event_p = dynamic_pointer_cast<Deletion>(events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime)));
			j_5_del_value_p = &(j_5_del_event_p->deletion_value);
		}
		else{j_5_del_value_p = &no_del_buffer;}

	}

}

void Coverage_err_counter::count_scenario(double scenario_seq_joint_proba , double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists){
	if(count_on_v){

		//Get mismatch list
		vector<int>& v_mismatch_list = *mismatches_lists.at(V_gene_seq);

		//Get the coverage
		//Get the length of the gene and a pointer to the right array to write on
		tmp_corr_len = v_gene_nucleotide_coverage_seq_p[**vgene_real_index_p].first;
		tmp_cov_p = v_gene_nucleotide_coverage_seq_p[**vgene_real_index_p].second;
		tmp_err_p = v_gene_per_nucleotide_error_seq_p[**vgene_real_index_p].second;

		/*
		 * Start at position 0
		 * Assume V is on the left of the read and compute left bound: max(0,-(**vgene_offset_p))
		 * Disregard P nucleotides, and set end bound as: tmp_corr_len - max(0,*v_3_del_value_p)
		 */

		this->recurs_coverage_count(scenario_seq_joint_proba,0,max(0,-(**vgene_offset_p)),tmp_corr_len - max(0,*v_3_del_value_p),tmp_corr_len);

		/*
		 * compute the position on the mismatch vector of the first Pnuc error and set it as the end bound
		 */
		size_t tmp_len_util = v_mismatch_list.size();
		for( i = 0 ; i != tmp_len_util ; ++i){
			//Disregard mismatches due to P nucleotides
			if(	 (v_mismatch_list[i]-(**vgene_offset_p))>=tmp_corr_len ){
				tmp_len_util = i;
				break;
			}
		}
		this->recurs_errors_count(scenario_seq_joint_proba,v_mismatch_list,vgene_offset_p,0,0,tmp_len_util,tmp_corr_len);

/*		//Get the corrected number of deletions(no negative deletion)
		tmp_corr_len -= max(0,*v_3_del_value_p); //FIXME assumes that V is on the left of the read

		// Compute the coverage
		for( i = max(0,-(**vgene_offset_p)) ; i != tmp_corr_len ; ++i ){
			tmp_cov_p[i]+=scenario_seq_joint_proba;
		}

		//Compute the error per nucleotide on the gene
		tmp_len_util = v_mismatch_list.size();
		for( i = 0 ; i != tmp_len_util ; ++i){
			//Disregard mismatches due to P nucleotides
			if(	 (v_mismatch_list[i]-(**vgene_offset_p))<tmp_corr_len ){
				tmp_err_p[v_mismatch_list[i]-(**vgene_offset_p)]+=scenario_seq_joint_proba;
			}
		}*/

	}

	if(count_on_d){
		throw("/!\\ D coverage and errors counters not implemented yet! /!\\ ");
	}

	if(count_on_j){

		//Get mismatch list
		vector<int>& j_mismatch_list = *mismatches_lists.at(J_gene_seq);

		//Get the coverage
		//Get the length of the gene and a pointer to the right array to write on
		tmp_corr_len = j_gene_nucleotide_coverage_seq_p[**jgene_real_index_p].first;
		tmp_cov_p = j_gene_nucleotide_coverage_seq_p[**jgene_real_index_p].second;
		tmp_err_p = j_gene_per_nucleotide_error_seq_p[**jgene_real_index_p].second;

		/*
		 * Start at position 0
		 * Assume V is on the left of the read and compute left bound: max(0,-(**vgene_offset_p))
		 * Disregard P nucleotides, and set end bound as: tmp_corr_len - max(0,*v_3_del_value_p)
		 */

		this->recurs_coverage_count(scenario_seq_joint_proba,0,max(0,(*j_5_del_value_p)),max(0,(*j_5_del_value_p))+(seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime) +1),tmp_corr_len);

		/*
		 * compute the position on the mismatch vector of the first Pnuc error and set it as the end bound
		 */
		size_t tmp_len_util = j_mismatch_list.size();
		for( i = 0 ; i != tmp_len_util ; ++i){
			//Disregard mismatches due to P nucleotides
			if(	 (j_mismatch_list[i] >= (**jgene_offset_p) + tmp_corr_len) ){
				tmp_len_util = i;
				break;
			}
		}
		this->recurs_errors_count(scenario_seq_joint_proba,j_mismatch_list,jgene_offset_p,0,tmp_len_util,j_mismatch_list.size(),tmp_corr_len);

/*		//Get the corrected number of deletions(no negative deletion)
		tmp_corr_len = max(0,(*j_5_del_value_p));

		// Compute the coverage
		const int tmp = (seq_offsets.at(J_gene_seq,Three_prime) - seq_offsets.at(J_gene_seq,Five_prime) +1);
		for( i = 0 ; i != tmp ; ++i ){
			tmp_cov_p[i+tmp_corr_len]+=scenario_seq_joint_proba;
		}

		//Compute the error per nucleotide on the gene
		tmp_len_util = j_mismatch_list.size();
		for( i = 0 ; i != tmp_len_util ; ++i){
			//Disregard mismatches due to P nucleotides
			if(	(j_mismatch_list[i] >= (**jgene_offset_p) + tmp_corr_len) ){
				tmp_err_p[j_mismatch_list[i]-(**jgene_offset_p)]+=scenario_seq_joint_proba;
			}
		}*/

	}
}

void Coverage_err_counter::count_sequence(double seq_likelihood , const Model_marginals& single_seq_marginals , const Model_Parms& single_seq_model_parms){
	//Normalize by the sequence likelihood and clean counter if not all seq are dumped
	if(seq_likelihood!=0){
		if(count_on_v){
			this->normalize_and_add_cov_and_err(seq_likelihood , n_v_real , v_gene_nucleotide_coverage_p , v_gene_per_nucleotide_error_p , v_gene_nucleotide_coverage_seq_p , v_gene_per_nucleotide_error_seq_p);
		}
		if(count_on_d){
			this->normalize_and_add_cov_and_err(seq_likelihood , n_d_real , d_gene_nucleotide_coverage_p , d_gene_per_nucleotide_error_p , d_gene_nucleotide_coverage_seq_p , d_gene_per_nucleotide_error_seq_p);
		}
		if(count_on_j){
			this->normalize_and_add_cov_and_err(seq_likelihood , n_j_real , j_gene_nucleotide_coverage_p , j_gene_per_nucleotide_error_p , j_gene_nucleotide_coverage_seq_p , j_gene_per_nucleotide_error_seq_p);
		}
	}
}

/*
 * Will copy multi sequence information on this
 * Will clean the other counter at the same time
 */
void Coverage_err_counter::add_checked(shared_ptr<Counter> counter){
	shared_ptr<Coverage_err_counter> other = dynamic_pointer_cast<Coverage_err_counter>(counter);

	//TODO add checks on counter nature and content
	double identity = 1.0;
	if(count_on_v){
		this->normalize_and_add_cov_and_err(identity , n_v_real , this->v_gene_nucleotide_coverage_p , this->v_gene_per_nucleotide_error_p , other->v_gene_nucleotide_coverage_p , other->v_gene_per_nucleotide_error_p);
	}
	if(count_on_d){
		this->normalize_and_add_cov_and_err(identity , n_d_real , this->d_gene_nucleotide_coverage_p , this->d_gene_per_nucleotide_error_p , other->d_gene_nucleotide_coverage_p , other->d_gene_per_nucleotide_error_p);
	}
	if(count_on_j){
		this->normalize_and_add_cov_and_err(identity , n_j_real , this->j_gene_nucleotide_coverage_p , this->j_gene_per_nucleotide_error_p , other->j_gene_nucleotide_coverage_p , other->j_gene_per_nucleotide_error_p);
	}
}

/*
 * Will output per sequence coverage and errors if needed
 * Also cleans individual seq counters at the same time
 */
void Coverage_err_counter::dump_sequence_data(int seq_index , int iteration_n){
	if(dump_individual_seqs){
		if(count_on_v){
			this->dump_cov_and_err_arrays(iteration_n,seq_index,output_cov_err_v_file_ptr,n_v_real,v_gene_nucleotide_coverage_seq_p,v_gene_per_nucleotide_error_seq_p);
		}

		if(count_on_d){
			this->dump_cov_and_err_arrays(iteration_n,seq_index,output_cov_err_d_file_ptr,n_d_real,d_gene_nucleotide_coverage_seq_p,d_gene_per_nucleotide_error_seq_p);
		}

		if(count_on_j){
			this->dump_cov_and_err_arrays(iteration_n,seq_index,output_cov_err_j_file_ptr,n_j_real,j_gene_nucleotide_coverage_seq_p,j_gene_per_nucleotide_error_seq_p);
		}
	}
}

void Coverage_err_counter::dump_data_summary(int iteration_n){
	if(not dump_individual_seqs){
		if(count_on_v){
			this->dump_cov_and_err_arrays(iteration_n,-1,output_cov_err_v_file_ptr,n_v_real,v_gene_nucleotide_coverage_p,v_gene_per_nucleotide_error_p);
		}

		if(count_on_d){
			this->dump_cov_and_err_arrays(iteration_n,-1,output_cov_err_d_file_ptr,n_d_real,d_gene_nucleotide_coverage_p,d_gene_per_nucleotide_error_p);
		}

		if(count_on_j){
			this->dump_cov_and_err_arrays(iteration_n,-1,output_cov_err_j_file_ptr,n_j_real,j_gene_nucleotide_coverage_p,j_gene_per_nucleotide_error_p);
		}
	}
}

shared_ptr<Counter> Coverage_err_counter::copy() const{
	shared_ptr<Coverage_err_counter> counter_copy_ptr (new Coverage_err_counter(path_to_file , count_on , record_Npoint_occurence , dump_individual_seqs , last_iter_only));
	counter_copy_ptr->fstreams_created = this->fstreams_created;
	if(this->fstreams_created){
		if(count_on_v){
			counter_copy_ptr->output_cov_err_v_file_ptr = this->output_cov_err_v_file_ptr;
		}
		if(count_on_d){
			counter_copy_ptr->output_cov_err_d_file_ptr = this->output_cov_err_d_file_ptr;
		}
		if(count_on_j){
			counter_copy_ptr->output_cov_err_j_file_ptr = this->output_cov_err_j_file_ptr;
		}
	}
	else{
		throw runtime_error("Counters should not be copied before stream initialization");
	}
	return counter_copy_ptr;
}

void Coverage_err_counter::allocate_coverage_and_errors_arrays(size_t n_real, const unordered_map<string , Event_realization> realizations ,pair<size_t,double*>*& gene_nucleotide_coverage_p,pair<size_t,double*>*& gene_per_nucleotide_error_p,pair<size_t,double*>*& gene_nucleotide_coverage_seq_p,pair<size_t,double*>*& gene_per_nucleotide_error_seq_p){

	//Create coverage and errors arrays
	gene_nucleotide_coverage_p = new pair<size_t,double*>[n_real];
	gene_per_nucleotide_error_p = new pair<size_t,double*>[n_real];
	gene_nucleotide_coverage_seq_p = new pair<size_t,double*>[n_real];
	gene_per_nucleotide_error_seq_p = new pair<size_t,double*>[n_real];

	for(unordered_map<string , Event_realization>::const_iterator iter = realizations.begin() ; iter != realizations.end() ; iter++){

		size_t tmp = pow((*iter).second.value_str_int.size(),record_Npoint_occurence);

		//Initialize normalized counters
		gene_nucleotide_coverage_p[(*iter).second.index] = pair<size_t,double*>((*iter).second.value_str_int.size(),new double [tmp]);
		gene_per_nucleotide_error_p[(*iter).second.index] = pair<size_t,double*>((*iter).second.value_str_int.size(),new double [tmp]);

		//Initialize sequence counters
		gene_nucleotide_coverage_seq_p[(*iter).second.index] = pair<size_t,double*>((*iter).second.value_str_int.size(),new double [tmp]);
		gene_per_nucleotide_error_seq_p[(*iter).second.index] = pair<size_t,double*>((*iter).second.value_str_int.size(),new double [tmp]);

		for(i=0 ; i!=pow((*iter).second.value_str_int.size(),record_Npoint_occurence) ; ++i){
			gene_nucleotide_coverage_p[(*iter).second.index].second[i]=0;
			gene_per_nucleotide_error_p[(*iter).second.index].second[i]=0;

			gene_nucleotide_coverage_seq_p[(*iter).second.index].second[i]=0;
			gene_per_nucleotide_error_seq_p[(*iter).second.index].second[i]=0;
		}
	}

}

void Coverage_err_counter::deallocate_coverage_and_errors_arrays(size_t n_real, const unordered_map<string , Event_realization> realizations ,pair<size_t,double*>*& gene_nucleotide_coverage_p,pair<size_t,double*>*& gene_per_nucleotide_error_p,pair<size_t,double*>*& gene_nucleotide_coverage_seq_p,pair<size_t,double*>*& gene_per_nucleotide_error_seq_p){

	if(n_real!=0){
		//If n_real==0 then the Counter has probably not been initialized
		for(unordered_map<string , Event_realization>::const_iterator iter = realizations.begin() ; iter != realizations.end() ; iter++){

			//Deallocate normalized counters
			delete [] gene_nucleotide_coverage_p[(*iter).second.index].second;
			delete [] gene_per_nucleotide_error_p[(*iter).second.index].second;

			//Deallocate sequence counters
			delete [] gene_nucleotide_coverage_seq_p[(*iter).second.index].second;
			delete [] gene_per_nucleotide_error_seq_p[(*iter).second.index].second;

		}

		//Deallocate coverage and errors arrays
		delete [] gene_nucleotide_coverage_p;
		delete [] gene_per_nucleotide_error_p;
		delete [] gene_nucleotide_coverage_seq_p;
		delete [] gene_per_nucleotide_error_seq_p;
	}
}

void Coverage_err_counter::dump_cov_and_err_arrays( int iteration_n ,  int seq_index , shared_ptr<ofstream> outfile_ptr , size_t n_real , pair<size_t,double*>* coverage_array_p , pair<size_t,double*>* error_array_p ){
	for(i=0 ; i!=n_real; ++i ){

		tmp_len_util = pow(coverage_array_p[i].first,record_Npoint_occurence);
		tmp_cov_p = coverage_array_p[i].second;
		tmp_err_p = error_array_p[i].second;

		if(dump_individual_seqs){
			(*outfile_ptr.get())<<iteration_n<<";"<<seq_index<<";"<<i<<";(";
		}
		else{
			(*outfile_ptr.get())<<iteration_n<<";"<<i<<";(";
		}

		//Symmetrize the array
		this->symmetrize_counter_array(tmp_cov_p,0,0,coverage_array_p[i].first);
		//Output it
		for(size_t j=0 ; j!=tmp_len_util ; ++j ){
			if(j!=0) (*outfile_ptr.get())<<",";
			(*outfile_ptr.get())<<tmp_cov_p[j];

			tmp_cov_p[j] = 0;
		}
		(*outfile_ptr.get())<<");(";

		//Symmetrize the array
		this->symmetrize_counter_array(tmp_err_p,0,0,coverage_array_p[i].first);
		//Output error array
		for(size_t j=0 ; j!=tmp_len_util ; ++j ){
			if(j!=0) (*outfile_ptr.get())<<",";
			(*outfile_ptr.get())<<tmp_err_p[j];

			tmp_err_p[j] = 0;
		}
		(*outfile_ptr.get())<<")"<<endl;
	}
}

void Coverage_err_counter::normalize_and_add_cov_and_err(double& normalizing_cst , size_t n_real , pair<size_t,double*>* target_coverage_array_p , pair<size_t,double*>* target_error_array_p , pair<size_t,double*>* base_coverage_array_p , pair<size_t,double*>* base_error_array_p){
	for(i=0 ; i!=n_real; ++i ){
		tmp_len_util = pow(base_coverage_array_p[i].first,record_Npoint_occurence);
		tmp_cov_p = base_coverage_array_p[i].second;
		tmp_err_p = base_error_array_p[i].second;

		for(size_t j=0 ; j!=tmp_len_util ; ++j){
			tmp_cov_p[j]/=normalizing_cst;
			tmp_err_p[j]/=normalizing_cst;
		}

		if(not dump_individual_seqs){
			for(size_t j=0 ; j!=tmp_len_util ; ++j){
				target_coverage_array_p[i].second[j]+=tmp_cov_p[j];
				tmp_cov_p[j] = 0; //Clean seq counter

				target_error_array_p[i].second[j]+=tmp_err_p[j];
				tmp_err_p[j] = 0; //Clean seq counter
			}
		}
	}
}

void Coverage_err_counter::recurs_coverage_count(double scenario_seq_joint_proba , size_t N , size_t begin_bound , size_t end_bound , size_t gene_len){
	for(size_t j = begin_bound ; j!=end_bound ; ++j){
		this->positions[N] = j;
		if(N<record_Npoint_occurence-1){
			this->recurs_coverage_count(scenario_seq_joint_proba , N+1 , j , end_bound , gene_len);
		}
		else{
			size_t adress = 0;
			for(size_t a = 0 ; a!=record_Npoint_occurence ; ++a){
				adress+=positions[a]*pow(gene_len,a);
			}
			tmp_cov_p[adress] += scenario_seq_joint_proba;
		}
	}
}

void Coverage_err_counter::recurs_errors_count(double scenario_seq_joint_proba , vector<int>& mismatch_list , const int** gene_offset_p , size_t N , size_t begin_bound , size_t end_bound , size_t gene_len){
	for(size_t j = begin_bound ; j!=end_bound ; ++j){
		this->positions[N] = mismatch_list.at(j)-(**gene_offset_p);
		if(N<record_Npoint_occurence-1){
			this->recurs_errors_count(scenario_seq_joint_proba , mismatch_list , gene_offset_p , N+1 , j , end_bound , gene_len);
		}
		else{
			size_t adress = 0;
			for(size_t a = 0 ; a!=record_Npoint_occurence ; ++a){
				adress+=positions[a]*pow(gene_len,a);
			}
			tmp_cov_p[adress] += scenario_seq_joint_proba;
		}
	}
}

void Coverage_err_counter::symmetrize_counter_array(double* counter_array , size_t N , size_t begin_bound , size_t gene_len){
	if(record_Npoint_occurence>1){
		for(size_t j = begin_bound ; j!=gene_len ; ++j){
			this->positions[N] = j;
			if(N<record_Npoint_occurence-1){
				this->symmetrize_counter_array(counter_array , N+1 , j  , gene_len);
			}
			else{
				size_t adress = 0;
				for(size_t a = 0 ; a!=record_Npoint_occurence ; ++a){
					adress+=positions[a]*pow(gene_len,a);
				}
				size_t* position_array = new size_t[record_Npoint_occurence];
				symmetrize_counter_array_recurs(adress , 0 , position_array , counter_array , gene_len);
				delete [] position_array;
			}
		}
	}
}

void Coverage_err_counter::symmetrize_counter_array_recurs(size_t adress , size_t N , size_t* position_array ,double* counter_array , size_t gene_len){
	for(size_t j=0 ; j!=record_Npoint_occurence ; ++j){

		bool is_valid = true;
		for(size_t k=0 ; k!= N ; ++k){
			if(j==position_array[k]){
				is_valid = false;
				break;
			}
		}

		if(is_valid){
			position_array[N] = j;
			if(N<record_Npoint_occurence-1){
				symmetrize_counter_array_recurs(adress,N+1,position_array,counter_array,gene_len);
			}
			else{
				size_t new_adress = 0;
				for(size_t a = 0 ; a!=record_Npoint_occurence ; ++a){
					new_adress+=positions[a]*pow(gene_len,position_array[a]);
				}
				counter_array[new_adress] = counter_array[adress];
			}
		}
	}
}
