/*
 * Deletion.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: quentin
 */

#include "Deletion.h"


using namespace std;
Deletion::Deletion(): Rec_Event() {
	this->type = Event_type::Deletion_t;
	this->update_event_name();
}

Deletion::Deletion(Gene_class gene, Seq_side del_side , pair<int,int> del_range):Rec_Event(gene , del_side , *(new unordered_map<string,Event_realization>)  ){
	//TODO prevent undefined_side entry(throw exception)

	int min_del = min(del_range.first , del_range.second);
	int max_del = max(del_range.first , del_range.second);

	this->type = Event_type::Deletion_t;
	this->len_max = -min_del;
	this->len_min = -max_del;

	for(int i=min_del ; i!=max_del+1 ; ++i){
		this->add_realization(i);
	}
	this->update_event_name();
}

Deletion::Deletion(Gene_class gene , Seq_side side , unordered_map<string,Event_realization>& realizations):Rec_Event(gene,side,realizations){
	this->type = Event_type::Deletion_t;
	for(unordered_map<string,Event_realization>::const_iterator iter = this->event_realizations.begin() ; iter != this->event_realizations.end() ; ++iter){
		if((*iter).second.value_int > (-this->len_min)){this->len_min = -(*iter).second.value_int;}
		else if((*iter).second.value_int < (-this->len_max)){this->len_max = -(*iter).second.value_int;}
	}
	this->update_event_name();
}

Deletion::~Deletion() {
	// TODO Auto-generated destructor stub
}

Rec_Event* Deletion::copy(){
	Deletion* new_deletion_p = new Deletion(this->event_class,this->event_side,this->event_realizations);//FIXME remove this new for all events and for the error rates
	new_deletion_p->priority = this->priority;
	new_deletion_p->nickname = this->nickname;
	new_deletion_p->fixed = this->is_fixed();
	new_deletion_p->update_event_name();
	new_deletion_p->set_event_identifier(this->event_index);
	return new_deletion_p;
}


void Deletion::add_realization(int del_number){
	this->Rec_Event::add_realization(*(new Event_realization(to_string(del_number) , del_number , "" , "" , this->size())));
	if(del_number > (-this->len_min)){this->len_min = (-del_number);}
	else if (del_number < (-this->len_max)){this->len_max = (-del_number);}
	this->update_event_name();
}
/**
 * General: Loop over all possible number of deletions for a given gene on a given sequence side
 *
 * Specific:
 * -First check whether any of these number of deletions is possible given the current position and number of deletions on other genes
 * -Loop over # of deletions in decreasing order
 */
void Deletion::iterate(double& scenario_proba , double& tmp_err_w_proba , const string& sequence , const string& int_sequence , Index_map& base_index_map , const unordered_map<Rec_Event_name,vector<pair<const Rec_Event*,int>>>& offset_map , queue<Rec_Event*>& model_queue , Marginal_array_p& updated_marginals_point , const Marginal_array_p& model_parameters_point ,const unordered_map<Gene_class , vector<Alignment_data>>& allowed_realizations , Seq_type_str_p_map& constructed_sequences , Seq_offsets_map& seq_offsets , Error_rate*& error_rate_p, const unordered_map<tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& events_map , Safety_bool_map& safety_set , Mismatch_vectors_map& mismatches_lists , double& seq_max_prob_scenario , double& proba_threshold_factor){

	base_index = base_index_map.at(this->event_index);
	//constructed_sequences_copy = constructed_sequences;
	//unordered_map<pair<Seq_type,Seq_side>,Seq_Offset> seq_offsets_copy (seq_offsets);
	//unordered_map<Rec_Event_name,int> base_index_map_copy(base_index_map);
	//unordered_map<Seq_type,vector<int>*> mismatches_lists_copy (mismatches_lists);



	switch((*this).event_class){

		case V_gene:
		{

			//v_3_offset = seq_offsets.at(pair<Seq_type,Seq_side>(V_gene_seq,Three_prime));
			v_3_offset = seq_offsets.at(V_gene_seq,Three_prime,memory_layer_offset_del-1);


			//Check D choice
			if(d_chosen){

				//if(safety_set.count(Event_safety::VD_safe) == 0){
				if(!safety_set.at(Event_safety::VD_safe,memory_layer_safety_1-1)){
					//d_5_offset = seq_offsets.at(pair<Seq_type,Seq_side>(D_gene_seq , Five_prime));
					//d_5_offset = seq_offsets.at(d_5_pair);
					d_5_offset = seq_offsets.at(D_gene_seq,Five_prime,memory_layer_offset_check1);

					d_5_min_offset = d_5_offset - d_5_min_del;
					d_5_max_offset = d_5_offset - d_5_max_del;

					vd_check = true;//Further check needed
				}
				else{
					vd_check = false;//No point of checking if D has not been picked because the offset is unknown
				}
			}
			else{
				vd_check = false;
			}



			//Check J choice
			if(j_chosen){

					//if(safety_set.count(Event_safety::VJ_safe) == 0){
				if(!safety_set.at(Event_safety::VJ_safe,memory_layer_safety_2-1)){
					//j_5_offset = seq_offsets.at(pair<Seq_type,Seq_side>(J_gene_seq , Five_prime));
					//j_5_offset = seq_offsets.at(j_5_pair);
					j_5_offset = seq_offsets.at(J_gene_seq,Five_prime,memory_layer_offset_check2);

					j_5_min_offset = j_5_offset - j_5_min_del;
					j_5_max_offset = j_5_offset - j_5_max_del;

					vj_check = true;//Further check needed
					}
				else{
					vj_check = false;//No point of checking if J has not been picked because the offset is unknown
				}
			}
			else{
				vj_check = false;
			}
			previous_str = (*constructed_sequences.at(V_gene_seq,memory_layer_cs-1));
			vector<int>& v_mismatch_list = *mismatches_lists.at(V_gene_seq , memory_layer_mismatches-1);


			for(forward_list<Event_realization>::const_iterator iter=(*this).int_value_and_index.begin() ; iter != (*this).int_value_and_index.end() ; ++iter){
				if((int)previous_str.size()>(*iter).value_int){ //Do not allow for deletion of the entire V
					//What about deletions going outside the read?
					//unordered_set<Event_safety> safety_set_copy = safety_set;

					v_3_new_offset = v_3_offset - (*iter).value_int;

					//There should be at least one nucleotide of the V in the read
					if(v_3_new_offset<0){
						continue;
					}

					if(vd_check){
						if( v_3_new_offset >= (d_5_max_offset)){
							//Even with maximum number of deletions on D overlap => bad alignments
							continue;
						}
						if( v_3_new_offset < (d_5_min_offset) ){
							//Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
							//safety_set_copy.emplace(Event_safety::VD_safe);
							safety_set.set_value(Event_safety::VD_safe,true,memory_layer_safety_1);
						}
						else{
							safety_set.set_value(Event_safety::VD_safe,false,memory_layer_safety_1);
						}
						//Already unsafe otherwise
					}

					if(vj_check){
						if( v_3_new_offset >= (j_5_max_offset)){
							//Even with maximum number of deletions on J overlap => bad alignments
							continue;
						}
						if( v_3_new_offset < (j_5_min_offset) ){
							//Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
							//safety_set_copy.emplace(Event_safety::VD_safe);
							//cout<<safety_set.get_current_memory_layer(Event_safety::VJ_safe)<<endl;
							safety_set.set_value(Event_safety::VJ_safe,true,memory_layer_safety_2);
						}
						else{
							safety_set.set_value(Event_safety::VJ_safe,false,memory_layer_safety_2);
						}
						//Already unsafe otherwise
					}




					new_index = base_index + (*iter).index;

					new_scenario_proba = scenario_proba;
					new_tmp_err_w_proba = tmp_err_w_proba;
					proba_contribution =1;

					iterate_common( iter , base_index_map , offset_map , model_parameters_point);


					//Positive or negative deletion (palindroms) mechanism
					if((*iter).value_int >= 0){
						//Delete the end of the V-gene (3' end)
						new_str = previous_str.substr(0,previous_str.size()-(*iter).value_int);

						//Update the mismatch list given the deletions
						if(!v_mismatch_list.empty()){
							mis_iter = v_mismatch_list.end();
							//Get an iterator referring to an actual integer (end() is not dereferenceable)
							mis_iter--;
							end_reached = false;
							while((*mis_iter)>v_3_new_offset ){
								if( mis_iter == v_mismatch_list.begin()){
									end_reached = true;
									break;
								}
								else{
									--mis_iter;
								}
							}

							if( end_reached){
								//Clear the vector but the capacity remains the same
								mismatches_vector.clear();
							}
							else{
								++mis_iter;
								mismatches_vector.assign(v_mismatch_list.begin(),mis_iter);
							}
						}
						else{
							mismatches_vector.clear();
						}
					}
					else //Negative deletions
					{
						if(v_3_new_offset < (int) sequence.size() ){
							//Check that the palindrom cannot be longer than the sequence itself
							if( (-(*iter).value_int) <= (int) previous_str.size()){
								//Copy the last nucleotides
								//cout<<"# p nucl: "<<(-(*iter).second.value_int)<<endl;
								tmp_str = previous_str.substr(previous_str.size()  + (*iter).value_int , string::npos );
								//cout<<"tmp_str: "<<tmp_str<<endl;
								//Reverse them
								//TODO revise this (perhaps not a good idea to perform these operations every time)
								reverse(tmp_str.begin() , tmp_str.end());
								make_transversions(tmp_str,true);
								//cout<<"rev_tmp_str: "<<tmp_str<<endl;
								//Merge strings
								new_str = previous_str + tmp_str;
								//cout<<"prev_str: "<<previous_str<<endl;
								//cout<<"new str:  "<<new_str<<endl;
								//Count mismatches and add them to the mismatches list
								mismatches_vector = v_mismatch_list;
							/*	cout<<"mismatch_vect : ";
								for(vector<int>::const_iterator test = mismatches_vector.begin() ; test != mismatches_vector.end() ; test++){
									cout<<(*test)<<";";
								}
								cout<<endl;
*/
								for(int i = 0 ; i != (-(*iter).value_int) ; ++i ){
									if(tmp_str[i] != int_sequence.at(v_3_offset+1+i)){
										mismatches_vector.push_back(v_3_offset+1+i);
									}
								}
/*								cout<<"prev_v_3_offset: "<<v_3_offset<<endl;
								cout<<"seq: "<<int_sequence<<endl;
								cout<<"new_mismatch_vect : ";
								for(vector<int>::const_iterator test = mismatches_vector.begin() ; test != mismatches_vector.end() ; test++){
									cout<<(*test)<<";";
								}
								cout<<endl;
								cout<<"--------------------------------------------------------------"<<endl;*/
							}
							else{
								continue;
							}
						}
						else{
							continue;
						}

					}

					constructed_sequences.set_value( V_gene_seq , &new_str , memory_layer_cs);
					//constructed_sequences_copy.at(V_gene_seq).erase(constructed_sequences.at(V_gene_seq).size() - (*iter).second.value_int);
					//Get rid of scenarios that delete more J nucleotides than the ones on the read //TODO improve this part (for J also)
					//if(constructed_sequences_copy.at(V_gene_seq).size()<1){continue;}//Already delt with upper

					//seq_offsets_copy.at(pair<Seq_type,Seq_side>(V_gene_seq,Three_prime)) = v_3_new_offset;
					//seq_offsets_copy.at(v_3_pair) = v_3_new_offset;
					seq_offsets.set_value(V_gene_seq,Three_prime,v_3_new_offset,memory_layer_offset_del);


					//Discard irrelevant mismatches (assuming the vector of mismatches is ordered) given the number of deletions

					mismatches_lists.set_value(V_gene_seq ,  &mismatches_vector , memory_layer_mismatches);

					//Add mismatches upper bound proba to the tmp_err_w_proba
					new_scenario_proba*=proba_contribution;
					new_tmp_err_w_proba*=pow(err_rate_upper_bound,mismatches_vector.size());
					compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
					if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
						//The order in which deletion are processed goes with decreasing number of deletion.
						//If a high number of deletions contains too many errors to be processed (even without taking the proba contribution into account), fewer deletions can only contain more thus the loop is broken
						break;
					}


					new_tmp_err_w_proba*=proba_contribution;
					compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
					if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
						continue;
					}


					Rec_Event::iterate_wrap_up(new_scenario_proba , new_tmp_err_w_proba , sequence , int_sequence , base_index_map , offset_map , model_queue  , updated_marginals_point , model_parameters_point , allowed_realizations , constructed_sequences , seq_offsets , error_rate_p , events_map , safety_set , mismatches_lists , seq_max_prob_scenario , proba_threshold_factor );
				}
			}
		}
			break;

		case D_gene:
			switch((*this).event_side){

			case Five_prime:
			{

				//d_5_offset = seq_offsets.at(pair<Seq_type,Seq_side>(D_gene_seq,Five_prime));
				//d_5_offset = seq_offsets.at(d_5_pair);
				d_5_offset = seq_offsets.at(D_gene_seq,Five_prime,memory_layer_offset_del-1);



				//Check V choice
				if(v_chosen){
					//if(safety_set.count(Event_safety::VD_safe) == 0){
					if(!safety_set.at(Event_safety::VD_safe,memory_layer_safety_1-1)){
						//v_3_offset = seq_offsets.at(pair<Seq_type,Seq_side>(V_gene_seq , Three_prime));
						//v_3_offset = seq_offsets.at(v_3_pair);

						v_3_offset = seq_offsets.at(V_gene_seq,Three_prime,memory_layer_offset_check1);

						v_3_max_offset = v_3_offset + v_3_min_del;
						v_3_min_offset = v_3_offset + v_3_max_del;

						vd_check = true;//Further check needed
					}
					else{
						vd_check = false;
					}

				}
				else{
					vd_check = false;
				}

				previous_str = (*constructed_sequences.at(D_gene_seq,memory_layer_cs-1));
				vector<int>& d_mismatch_list = *mismatches_lists.at(D_gene_seq , memory_layer_mismatches-1);


				for(forward_list<Event_realization>::const_iterator iter=(*this).int_value_and_index.begin() ; iter != (*this).int_value_and_index.end() ; ++iter){
					if( (int) previous_str.size()>=(*iter).value_int){

						//unordered_set<Event_safety> safety_set_copy = safety_set;

						d_5_new_offset = d_5_offset+(*iter).value_int;
						if(vd_check){

							if( d_5_new_offset <= (v_3_min_offset)){
								//Even with maximum number of deletions on V , D overlap => bad alignments
								continue;
							}
							if( d_5_new_offset > (v_3_max_offset) ){
								//Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
								//safety_set_copy.emplace(Event_safety::VD_safe);
								safety_set.set_value(Event_safety::VD_safe,true,memory_layer_safety_1);
							}
							else{
								safety_set.set_value(Event_safety::VD_safe,false,memory_layer_safety_1);
							}
							//Already unsafe otherwise
						}



						new_index = base_index + (*iter).index;
						new_scenario_proba = scenario_proba;
						new_tmp_err_w_proba = tmp_err_w_proba;
						proba_contribution = 1;

						iterate_common( iter , base_index_map , offset_map , model_parameters_point);

						//Positive or negative deletion (palindroms) mechanism
						if((*iter).value_int >= 0){

							//Delete the beginning of the D-gene (5' end)
							new_str = previous_str.substr((*iter).value_int,string::npos);

							//Update the mismatches given the number of deleted nucleotides
							//Discard irrelevant mismatches (assuming the vector of mismatches is ordered) given the number of deletions
							if(!d_mismatch_list.empty()){
								mis_iter = d_mismatch_list.begin();

								end_reached = false;
								while( ( (*mis_iter)< d_5_new_offset ) ){
									++mis_iter;
									if(mis_iter==d_mismatch_list.end()){
										end_reached = true;
										break;
									}
								}
								if(end_reached){
									mismatches_vector.clear();
								}
								else{
									mismatches_vector.assign(mis_iter,d_mismatch_list.end());
								}
							}
							else{
								mismatches_vector.clear();
							}

						}else //Negative deletion
						{
							if(d_5_new_offset>=0){
								//Check that the palindrom cannot be longer than the sequence itself
								if( (-(*iter).value_int) <= (int) previous_str.size()){

									//Copy the first nucleotides
									tmp_str = previous_str.substr(0 , -(*iter).value_int );

									//Reverse them
									reverse(tmp_str.begin() , tmp_str.end());
									make_transversions(tmp_str,true);
									//Merge strings
									new_str =  tmp_str + previous_str ;

									mismatches_vector = d_mismatch_list;

									for(int i = 0 ; i != (-(*iter).value_int) ; ++i ){
										//Check for errors in reverse order to keep the mismatch vector ordered
										if(tmp_str[i] != int_sequence.at(d_5_new_offset+i)){
											mismatches_vector.push_back(d_5_new_offset+i);
										}
									}
									sort(mismatches_vector.begin(),mismatches_vector.end());

								}
								else{
									continue;
								}
							}
							else{
								continue;
							}

						}


						constructed_sequences.set_value(D_gene_seq , &new_str , memory_layer_cs);
						//constructed_sequences_copy.at(D_gene_seq).erase(0 , (*iter).second.value_int);

						//seq_offsets_copy.at(pair<Seq_type,Seq_side>(D_gene_seq,Five_prime)) = d_5_new_offset;
						//seq_offsets_copy.at(d_5_pair) = d_5_new_offset;
						seq_offsets.set_value(D_gene_seq,Five_prime,d_5_new_offset,memory_layer_offset_del);




						mismatches_lists.set_value(D_gene_seq , &mismatches_vector , memory_layer_mismatches);
						//TODO add mismatches if del_d3 has been processed
						new_scenario_proba*=proba_contribution;
						if(d_del_opposite_side_processed){
							new_tmp_err_w_proba*=pow(err_rate_upper_bound,mismatches_vector.size());
							compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
							if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
								//The order in which deletion are processed goes with decreasing number of deletion.
								//If a high number of deletions contains too many errors to be processed (even without taking the proba contribution into account), fewer deletions can only contain more thus the loop is broken
								break;
							}


							new_tmp_err_w_proba*=proba_contribution;
							compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
							if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
								continue;
							}
						}
						else{
							new_tmp_err_w_proba*=proba_contribution;
							compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
							if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
								continue;
							}
						}



						Rec_Event::iterate_wrap_up(new_scenario_proba , new_tmp_err_w_proba , sequence , int_sequence , base_index_map , offset_map , model_queue  , updated_marginals_point  , model_parameters_point , allowed_realizations , constructed_sequences , seq_offsets , error_rate_p , events_map , safety_set , mismatches_lists ,seq_max_prob_scenario , proba_threshold_factor);

					}
				}
			}

				break;

			case Three_prime:
			{


				//d_3_offset = seq_offsets.at(pair<Seq_type,Seq_side>(D_gene_seq,Three_prime));
				//d_3_offset = seq_offsets.at(d_3_pair);
				d_3_offset = seq_offsets.at(D_gene_seq,Three_prime,memory_layer_offset_del-1);


				//Check J choice
				if(j_chosen){
					//if(safety_set.count(Event_safety::DJ_safe) == 0){
					if(!safety_set.at(Event_safety::DJ_safe,memory_layer_safety_2-1)){
						//j_5_offset = seq_offsets.at(pair<Seq_type,Seq_side>(J_gene_seq , Five_prime));
						//j_5_offset = seq_offsets.at(j_5_pair);
						j_5_offset = seq_offsets.at(J_gene_seq,Five_prime,memory_layer_offset_check2);

						j_5_min_offset = j_5_offset - j_5_min_del;
						j_5_max_offset = j_5_offset - j_5_max_del;

						dj_check = true;//Further check needed
					}
					else{
						dj_check = false;
					}
				}
				else{
					dj_check = false;
				}

				previous_str = (*constructed_sequences.at(D_gene_seq,memory_layer_cs-1));
				vector<int>& d_mismatch_list = *mismatches_lists.at(D_gene_seq , memory_layer_mismatches-1);

				for(forward_list<Event_realization>::const_iterator iter=(*this).int_value_and_index.begin() ; iter != (*this).int_value_and_index.end() ; ++iter){
					if((int)previous_str.size()>=(*iter).value_int){

						//unordered_set<Event_safety> safety_set_copy = safety_set;

						d_3_new_offset = d_3_offset-(*iter).value_int;
						if(dj_check){
							if( d_3_new_offset >= (j_5_max_offset) ){
								//Even with maximum number of deletions on J overlap => bad alignments
								continue;
							}
							if( d_3_new_offset < (j_5_min_offset) ){
								//Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
								//safety_set_copy.emplace(Event_safety::VJ_safe);
								safety_set.set_value(Event_safety::DJ_safe,true,memory_layer_safety_2);
							}
							else{
								safety_set.set_value(Event_safety::DJ_safe,false,memory_layer_safety_2);
							}
							//Already unsafe otherwise
						}


						new_index = base_index + (*iter).index;
						new_scenario_proba = scenario_proba;
						new_tmp_err_w_proba = tmp_err_w_proba;
						proba_contribution = 1;

						iterate_common( iter , base_index_map , offset_map , model_parameters_point);


						//Positive or negative deletion (palindroms) mechanism
						if((*iter).value_int >= 0){
							//Delete the end of the D-gene (3'end)
							new_str = previous_str.substr(0,previous_str.size() - (*iter).value_int);

							//Update mismatches list given the number of deletions
							//Discard irrelevant mismatches (assuming the vector of mismatches is ordered) given the number of deletions
							if(!d_mismatch_list.empty()){
								mis_iter = d_mismatch_list.end();
								//Get an iterator referring to an actual integer (end() is not dereferenceable)
								mis_iter--;
								end_reached = false;
								while((*mis_iter)>d_3_new_offset ){
									if( mis_iter == d_mismatch_list.begin()){
										end_reached = true;
										break;
									}
									else{
										mis_iter--;
									}
								}

								if( end_reached){
									//Clear the vector but the capacity remains the same
									mismatches_vector.clear();
								}
								else{
									++mis_iter;
									mismatches_vector.assign(d_mismatch_list.begin(),mis_iter);
								}
							}
							else{
								mismatches_vector.clear();
							}

						}
						else //Negative deletion
						{
							if(d_3_new_offset < (int) sequence.size() ){
								//Check that the palindrom cannot be longer than the sequence itself
								if( (-(*iter).value_int) <= (int) previous_str.size()){
									//Copy the last nucleotides

									tmp_str = previous_str.substr(previous_str.size()  + (*iter).value_int , string::npos );

									//Reverse them
									reverse(tmp_str.begin() , tmp_str.end());
									make_transversions(tmp_str,true);

									//Merge strings
									new_str = previous_str + tmp_str;

									//Count mismatches and add them to the mismatches list
									mismatches_vector = d_mismatch_list;

									for(int i = 0 ; i != (-(*iter).value_int) ; ++i ){
										if(d_3_offset+1+i<0){
											cout<<"problem"<<endl;
											cout<<d_3_offset<<endl;
											cout<<i<<endl;

											cout<<seq_offsets.at(V_gene_seq,Three_prime)<<endl;
											cout<<seq_offsets.at(D_gene_seq,Five_prime)<<endl;
											cout<<seq_offsets.at(D_gene_seq,Three_prime)<<endl;
											cout<<seq_offsets.at(J_gene_seq,Five_prime)<<endl;

										}

										if(tmp_str[i] != int_sequence.at(d_3_offset+1+i)){
											mismatches_vector.push_back(d_3_offset+1+i);
										}
									}
								}
								else{
									continue;
								}
							}
							else{
								continue;
							}

						}


						constructed_sequences.set_value(D_gene_seq , &new_str,memory_layer_cs);
						//constructed_sequences_copy.at(D_gene_seq).erase(constructed_sequences.at(D_gene_seq).size() - (*iter).second.value_int);

						//seq_offsets_copy.at(pair<Seq_type,Seq_side>(D_gene_seq,Three_prime)) = d_3_new_offset;
						//seq_offsets_copy.at(d_3_pair) = d_3_new_offset;
						seq_offsets.set_value(D_gene_seq,Three_prime,d_3_new_offset,memory_layer_offset_del);






						mismatches_lists.set_value(D_gene_seq , &mismatches_vector , memory_layer_mismatches);
						//TODO add mismatches if ddel5 processed
						new_scenario_proba*=proba_contribution;
						if(d_del_opposite_side_processed){
							new_tmp_err_w_proba*=pow(err_rate_upper_bound,mismatches_vector.size());
							compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
							if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
								//The order in which deletion are processed goes with decreasing number of deletion.
								//If a high number of deletions contains too many errors to be processed (even without taking the proba contribution into account), fewer deletions can only contain more thus the loop is broken
								break;
							}


							new_tmp_err_w_proba*=proba_contribution;
							compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
							if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
								continue;
							}
						}
						else{
							new_tmp_err_w_proba*=proba_contribution;
							compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
							if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
								continue;
							}
						}

						Rec_Event::iterate_wrap_up(new_scenario_proba , new_tmp_err_w_proba , sequence , int_sequence , base_index_map , offset_map , model_queue  , updated_marginals_point  , model_parameters_point , allowed_realizations , constructed_sequences , seq_offsets , error_rate_p , events_map , safety_set , mismatches_lists ,seq_max_prob_scenario , proba_threshold_factor);

					}


				}
			}
				break;

			default:
				throw invalid_argument("Unknown side for D deletion: " + (*this).event_side); //TODO explicitly throw the side used
				break;

			}

			break;

		case J_gene:
		{

			//j_5_offset =  seq_offsets.at(pair<Seq_type,Seq_side>(J_gene_seq,Five_prime));
			//j_5_offset =  seq_offsets.at(j_5_pair);
			j_5_offset =  seq_offsets.at(J_gene_seq,Five_prime,memory_layer_offset_del-1);


			//Check D choice
			if(d_chosen){
				//if(safety_set.count(Event_safety::DJ_safe) == 0){
				if(!safety_set.at(Event_safety::DJ_safe,memory_layer_safety_2-1)){
					//d_3_offset = seq_offsets.at(pair<Seq_type,Seq_side>(D_gene_seq , Three_prime));
					//d_3_offset = seq_offsets.at(d_3_pair);
					d_3_offset = seq_offsets.at(D_gene_seq,Three_prime,memory_layer_offset_check2);

					d_3_min_offset = d_3_offset + d_3_max_del;
					d_3_max_offset = d_3_offset + d_3_min_del;

					dj_check = true;//Further check needed
				}
				else{
					dj_check = false;
				}
			}
			else{
				dj_check = false;
			}



			//Check V choice
			if(v_chosen){
				//if(safety_set.count(Event_safety::VJ_safe) == 0){
				if(!safety_set.at(Event_safety::VJ_safe,memory_layer_safety_1-1)){
					//v_3_offset = seq_offsets.at(pair<Seq_type,Seq_side>(V_gene_seq , Three_prime));
					//v_3_offset = seq_offsets.at(v_3_pair);
					v_3_offset = seq_offsets.at(V_gene_seq,Three_prime,memory_layer_offset_check1);

					v_3_min_offset = v_3_offset + v_3_max_del;
					v_3_max_offset = v_3_offset + v_3_min_del;

					vj_check = true;//Further check needed
				}
				else{
					vj_check = false;
				}
			}
			else{
				vj_check = false;
			}

			previous_str = (*constructed_sequences.at(J_gene_seq,memory_layer_cs-1));
			vector<int>& j_mismatch_list = *mismatches_lists.at(J_gene_seq , memory_layer_mismatches-1);
			for(forward_list<Event_realization>::const_iterator iter=(*this).int_value_and_index.begin() ; iter != (*this).int_value_and_index.end() ; ++iter){
				if( (int) previous_str.size()>(*iter).value_int){

					//unordered_set<Event_safety> safety_set_copy = safety_set;

					j_5_new_offset = j_5_offset +(*iter).value_int;
					if(vj_check){
						if( j_5_new_offset <= (v_3_min_offset)){
							//Even with maximum number of deletions on V overlap => bad alignments
							continue;
						}
						if( j_5_new_offset > (v_3_max_offset) ){
							//Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
							//safety_set_copy.emplace(Event_safety::VD_safe);
							safety_set.set_value(Event_safety::VJ_safe,true,memory_layer_safety_1);
						}
						else{
							safety_set.set_value(Event_safety::VJ_safe,false,memory_layer_safety_1);
						}
						//Already unsafe otherwise
					}
					if(dj_check){
						if( j_5_new_offset <= (d_3_min_offset) ){
							//Even with maximum number of deletions on each side the D and J overlap => bad alignments
							continue;
						}
						if( j_5_new_offset > (d_3_max_offset) ){
							//Even with minimum number of deletions there's no overlap => safe even without knowing the number of deletions
							//safety_set_copy.emplace(Event_safety::VJ_safe);
							safety_set.set_value(Event_safety::DJ_safe,true,memory_layer_safety_2);
						}
						else{
							safety_set.set_value(Event_safety::DJ_safe,false,memory_layer_safety_2);
						}
						//Already unsafe otherwise
					}



					new_index = base_index + (*iter).index;
					new_scenario_proba=scenario_proba;
					new_tmp_err_w_proba=tmp_err_w_proba;
					proba_contribution = 1;

					iterate_common( iter , base_index_map , offset_map , model_parameters_point);

					//Positive or negative deletion (palindroms) mechanism
					if((*iter).value_int >= 0){
						//Delete the beginning of the J-gene (5' end)
						new_str = previous_str.substr((*iter).value_int);

						//Update mismatch list given the number of deletions
						//Discard irrelevant mismatches (assuming the vector of mismatches is ordered) given the number of deletions
						if(!j_mismatch_list.empty()){
							mis_iter = j_mismatch_list.begin();
							end_reached = false;
							while( ( (*mis_iter)< j_5_new_offset ) ){
								++mis_iter;
								if(mis_iter==j_mismatch_list.end()){
									end_reached = true;
									break;
								}
							}
							if(end_reached){
								mismatches_vector.clear();
							}
							else{
								mismatches_vector.assign(mis_iter,j_mismatch_list.end());
							}
						}
						else{
							mismatches_vector.clear();
						}


					}
					else //Negative deletion
					{
						//Check that the palindrom cannot be longer than the sequence itself
						if( (-(*iter).value_int) <= (int)previous_str.size()){

							//Copy the first nucleotides
							//cout<<"# p nucl: "<<(-(*iter).second.value_int)<<endl;
							tmp_str = previous_str.substr(0 , -(*iter).value_int );
							//cout<<"tmp_str: "<<tmp_str<<endl;
							//Reverse them
							reverse(tmp_str.begin() , tmp_str.end());
							make_transversions(tmp_str,true);
							//Merge strings
							new_str =  tmp_str + previous_str ;

							//cout<<"prev_str: "<<previous_str<<endl;
							//cout<<"new str:  "<<new_str<<endl;
							//Count mismatches and add them to the mismatches list
							mismatches_vector = j_mismatch_list;
/*							cout<<"mismatch_vect : ";
							for(vector<int>::const_iterator test = mismatches_vector.begin() ; test != mismatches_vector.end() ; test++){
								cout<<(*test)<<";";
							}
							cout<<endl;*/
							for(int i = 0 ; i != (-(*iter).value_int) ; ++i ){
								//Check for errors in reverse order to keep the mismatch vector ordered
								if(tmp_str[i] != int_sequence.at(j_5_new_offset+i)){
									mismatches_vector.push_back(j_5_new_offset+i);
								}
							}
							sort(mismatches_vector.begin(),mismatches_vector.end());
/*							cout<<"prev_v_3_offset: "<<v_3_offset<<endl;
							cout<<"seq: "<<int_sequence<<endl;
							cout<<"new_mismatch_vect : ";
							for(vector<int>::const_iterator test = mismatches_vector.begin() ; test != mismatches_vector.end() ; test++){
								cout<<(*test)<<";";
							}
							cout<<endl;
							cout<<endl;*/
						}
						else{
							continue;
						}
					}



					constructed_sequences.set_value(J_gene_seq , &new_str , memory_layer_cs);
					//constructed_sequences_copy.at(J_gene_seq).erase(0 , (*iter).second.value_int);
					//Get rid of scenarios that delete more J nucleotides than the ones on the read //TODO improve this part (for J also)
					//if(constructed_sequences_copy.at(J_gene_seq).size()<1){continue;}

					//seq_offsets_copy.at(pair<Seq_type,Seq_side>(J_gene_seq,Five_prime)) = j_5_new_offset;
					//seq_offsets_copy.at(j_5_pair) = j_5_new_offset;
					seq_offsets.set_value(J_gene_seq,Five_prime, j_5_new_offset, memory_layer_offset_del);




					mismatches_lists.set_value(J_gene_seq , &mismatches_vector , memory_layer_mismatches);

					//Add mismatches upper bound proba to the tmp_err_w_proba
					new_scenario_proba*=proba_contribution;
					new_tmp_err_w_proba*=pow(err_rate_upper_bound,mismatches_vector.size());
					compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
					if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
						//The order in which deletion are processed goes with decreasing number of deletion.
						//If a high number of deletions contains too many errors to be processed (even without taking the proba contribution into account), fewer deletions can only contain more thus the loop is broken
						break;
					}


					new_tmp_err_w_proba*=proba_contribution;
					compute_upper_bound_scenario_proba(new_tmp_err_w_proba);
					if(scenario_upper_bound_proba<(seq_max_prob_scenario*proba_threshold_factor)){
						continue;
					}

					Rec_Event::iterate_wrap_up(new_scenario_proba , new_tmp_err_w_proba , sequence , int_sequence , base_index_map , offset_map , model_queue  , updated_marginals_point  , model_parameters_point , allowed_realizations , constructed_sequences , seq_offsets , error_rate_p , events_map , safety_set , mismatches_lists , seq_max_prob_scenario , proba_threshold_factor);
				}
			}
		}

			break;

		default:
			throw invalid_argument("Unknown gene for deletions : " + this->event_class);
			break;
	}

}


/*
 * This short method performs the iterate operations common to all Rec_event (modify index map and fetch realization probability)
 */
 void Deletion::iterate_common(forward_list<Event_realization>::const_iterator& iter , Index_map& base_index_map , const unordered_map<Rec_Event_name,vector<pair<const Rec_Event*,int>>>& offset_map ,const Marginal_array_p& model_parameters_point){
	 //realization_event_index  =

	 	/*if (offset_map.count(this->name)!=0){
	 		for(vector<pair<const Rec_Event*,int>>::const_iterator jiter = offset_map.at(this->name).begin() ; jiter!= offset_map.at(this->name).end() ; jiter++){
				//modify index map using offset map
				base_index_map_copy.at((*jiter).first->get_name()) = base_index_map.at((*jiter).first->get_name()) + (*iter).second.index*(*jiter).second;
	 		}
	 	}*/

	 for(forward_list<tuple<int,int,int>>::const_iterator jiter = memory_and_offsets.begin() ; jiter!=memory_and_offsets.end() ; ++jiter){
		//Get previous index for the considered event
		previous_marginal_index = base_index_map.at(get<0>(*jiter),get<1>(*jiter)-1);
		//Update the index given the realization and the offset
		previous_marginal_index += (*iter).index*get<2>(*jiter);
		//Set the value
		base_index_map.set_value(get<0>(*jiter) , previous_marginal_index , get<1>(*jiter));
	}


	 //Compute the probability of the scenario considering the realization (*iter) we're looking at
	 proba_contribution = (model_parameters_point[base_index+(*iter).index]);
 }

 queue<int> Deletion::draw_random_realization( const Marginal_array_p model_marginals_p , unordered_map<Rec_Event_name,int>& index_map , const unordered_map<Rec_Event_name,vector<pair<const Rec_Event*,int>>>& offset_map , unordered_map<Seq_type , string>& constructed_sequences , default_random_engine& generator)const{

	 uniform_real_distribution<double> distribution(0.0,1.0);
	 double rand = distribution(generator);
	double prob_count = 0;
	queue<int> realization_queue ;
	for(unordered_map<string,Event_realization>::const_iterator iter = this->event_realizations.begin() ; iter != this->event_realizations.end() ; ++iter ){
		prob_count += model_marginals_p[index_map.at(this->get_name()) + (*iter).second.index];
		if(prob_count>=rand){
			switch(this->event_class){

			case V_gene:
				if((*iter).second.value_int>=0){
					constructed_sequences.at(V_gene_seq).erase(constructed_sequences.at(V_gene_seq).size() - (*iter).second.value_int);
				}
				else{
					string& v_gene_seq = constructed_sequences.at(V_gene_seq);
					tmp_str = v_gene_seq.substr(v_gene_seq.size() + (*iter).second.value_int , string::npos);
					reverse(tmp_str.begin(),tmp_str.end());
					make_transversions(tmp_str,false);
					v_gene_seq+=tmp_str;
				}

				break;

			case D_gene:
				switch(this->event_side){

				case Five_prime:
					if((*iter).second.value_int>=0){
						constructed_sequences.at(D_gene_seq).erase(0 , (*iter).second.value_int);
					}
					else{
						string& d_gene_seq = constructed_sequences.at(D_gene_seq);
						tmp_str = d_gene_seq.substr(0 , -(*iter).second.value_int );
						reverse(tmp_str.begin(),tmp_str.end());
						make_transversions(tmp_str,false);
						new_str = tmp_str + d_gene_seq;
						d_gene_seq = new_str;

					}

					break;

				case Three_prime:
					if((*iter).second.value_int>=0){
						constructed_sequences.at(D_gene_seq).erase(constructed_sequences.at(D_gene_seq).size() - (*iter).second.value_int);
					}
					else{
						string& d_gene_seq = constructed_sequences.at(D_gene_seq);
						tmp_str = d_gene_seq.substr(d_gene_seq.size() + (*iter).second.value_int , string::npos);
						reverse(tmp_str.begin(),tmp_str.end());
						make_transversions(tmp_str,false);
						d_gene_seq+=tmp_str;

					}

					break;

				default:
					break;
				}
				break;
			case J_gene:
				if((*iter).second.value_int>=0){
					constructed_sequences.at(J_gene_seq).erase(0 , (*iter).second.value_int);
				}
				else{
					string& j_gene_seq = constructed_sequences.at(J_gene_seq);
					tmp_str = j_gene_seq.substr(0 , -(*iter).second.value_int );
					reverse(tmp_str.begin(),tmp_str.end());
					make_transversions(tmp_str,false);
					new_str = tmp_str + j_gene_seq;
					j_gene_seq = new_str;
				}

				break;
			default:
				break;

			}
			realization_queue.push((*iter).second.index);
			if(offset_map.count(this->get_name()) != 0){
				for (vector<pair<const Rec_Event*,int>>::const_iterator jiter = offset_map.at(this->get_name()).begin() ; jiter!= offset_map.at(this->get_name()).end() ; ++jiter){
					index_map.at((*jiter).first->get_name()) += (*iter).second.index*(*jiter).second;
				}
			}

			break;
		}
	}
	return realization_queue;
}
 void Deletion::write2txt(ofstream& outfile){
 	outfile<<"#Deletion;"<<event_class<<";"<<event_side<<";"<<priority<<";"<<nickname<<endl;
 	for(unordered_map<string,Event_realization>::const_iterator iter=event_realizations.begin() ; iter!= event_realizations.end() ; ++iter){
 		outfile<<"%"<<(*iter).second.value_int<<";"<<(*iter).second.index<<endl;
 	}
 }

 void Deletion::initialize_event( unordered_set<Rec_Event_name>& processed_events , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, Rec_Event*>& events_map , const unordered_map<Rec_Event_name,vector<pair<const Rec_Event*,int>>>& offset_map , Seq_type_str_p_map& constructed_sequences , Safety_bool_map& safety_set , Error_rate* error_rate_p , Mismatch_vectors_map& mismatches_list , Seq_offsets_map& seq_offsets , Index_map& index_map){


	 err_rate_upper_bound = error_rate_p->get_err_rate_upper_bound();


	 //TODO change this and the usage of int_value_and_index
	 int_value_and_index.clear();
	 for(unordered_map<string,Event_realization>::const_iterator iter=(*this).event_realizations.begin() ; iter != (*this).event_realizations.end() ; ++iter){
		 int_value_and_index.push_front((*iter).second);
	 }
	 int_value_and_index.sort(del_numb_compare);


	 //Check V choice
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side))!=0){
		const Rec_Event* v_choice_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,V_gene,Undefined_side));
		if(processed_events.count(v_choice_p->get_name())!=0){v_chosen = true;}
		else{v_chosen=false;}
	}
	else{v_chosen=false;}

	//Check D choice
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side))!=0){
		const Rec_Event* d_choice_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,D_gene,Undefined_side));
		if(processed_events.count(d_choice_p->get_name())!=0){d_chosen = true;}
		else{d_chosen=false;}
	}
	else{d_chosen=false;}

	//Check J choice
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side))!=0){
		const Rec_Event* j_choice_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(GeneChoice_t,J_gene,Undefined_side));
		if(processed_events.count(j_choice_p->get_name())!=0){j_chosen = true;}
		else{j_chosen=false;}
	}
	else{j_chosen=false;}

	switch(this->event_class){
			case V_gene:
				seq_offsets.request_memory_layer(V_gene_seq,Three_prime);
				memory_layer_offset_del = seq_offsets.get_current_memory_layer(V_gene_seq,Three_prime);
				mismatches_list.request_memory_layer(V_gene_seq);
				this->memory_layer_mismatches = mismatches_list.get_current_memory_layer(V_gene_seq);
				constructed_sequences.request_memory_layer(V_gene_seq);
				this->memory_layer_cs = constructed_sequences.get_current_memory_layer(V_gene_seq);
				if(d_chosen){
					safety_set.request_memory_layer(VD_safe);
					memory_layer_safety_1 = safety_set.get_current_memory_layer(VD_safe);
					memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(D_gene_seq,Five_prime);
					//cout<<"v_del_1 : "<<memory_layer_safety_1<<endl;
				}
				if(j_chosen){
					safety_set.request_memory_layer(VJ_safe);
					memory_layer_safety_2 = safety_set.get_current_memory_layer(VJ_safe);
					memory_layer_offset_check2 = seq_offsets.get_current_memory_layer(J_gene_seq,Five_prime);
					//cout<<"v_del_2 : "<<memory_layer_safety_2<<endl;
				}

				break;
			case D_gene:
				mismatches_list.request_memory_layer(D_gene_seq);
				this->memory_layer_mismatches = mismatches_list.get_current_memory_layer(D_gene_seq);
				constructed_sequences.request_memory_layer(D_gene_seq);
				this->memory_layer_cs = constructed_sequences.get_current_memory_layer(D_gene_seq);
				switch(this->event_side){
					case Five_prime:
						seq_offsets.request_memory_layer(D_gene_seq,Five_prime);
						memory_layer_offset_del = seq_offsets.get_current_memory_layer(D_gene_seq,Five_prime);
						if(v_chosen){
							safety_set.request_memory_layer(VD_safe);
							memory_layer_safety_1 = safety_set.get_current_memory_layer(VD_safe);
							memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(V_gene_seq,Three_prime);
							//cout<<"d_del_1: "<<memory_layer_safety_1<<endl;
						}
						if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime)) != 0){
							const Rec_Event* del_d_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime));
							if(processed_events.count(del_d_p->get_name())!=0){
								d_del_opposite_side_processed = true;
							}
							else{
								d_del_opposite_side_processed = false;
							}
						}
						else{
							d_del_opposite_side_processed = true;
						}
						break;
					case Three_prime:
						seq_offsets.request_memory_layer(D_gene_seq,Three_prime);
						memory_layer_offset_del = seq_offsets.get_current_memory_layer(D_gene_seq,Three_prime);
						if(j_chosen){
							safety_set.request_memory_layer(DJ_safe);
							memory_layer_safety_2 = safety_set.get_current_memory_layer(DJ_safe);
							memory_layer_offset_check2 = seq_offsets.get_current_memory_layer(J_gene_seq,Five_prime);
							//cout<<"d_del_2: "<<memory_layer_safety_2<<endl;
						}
						if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime)) != 0){
							const Rec_Event* del_d_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime));
							if(processed_events.count(del_d_p->get_name())!=0){
								d_del_opposite_side_processed = true;
							}
							else{
								d_del_opposite_side_processed = false;
							}
						}
						else{
							d_del_opposite_side_processed = true;
						}
						break;
				}

				break;
			case J_gene:
				seq_offsets.request_memory_layer(J_gene_seq,Five_prime);
				memory_layer_offset_del = seq_offsets.get_current_memory_layer(J_gene_seq,Five_prime);
				mismatches_list.request_memory_layer(J_gene_seq);
				this->memory_layer_mismatches = mismatches_list.get_current_memory_layer(J_gene_seq);
				constructed_sequences.request_memory_layer(J_gene_seq);
				this->memory_layer_cs = constructed_sequences.get_current_memory_layer(J_gene_seq);
				if(v_chosen){
					safety_set.request_memory_layer(VJ_safe);
					memory_layer_safety_1 = safety_set.get_current_memory_layer(VJ_safe);
					memory_layer_offset_check1 = seq_offsets.get_current_memory_layer(V_gene_seq,Three_prime);
					//cout<<"j_del_1: "<<memory_layer_safety_1<<endl;
				}
				if(d_chosen){
					safety_set.request_memory_layer(DJ_safe);
					memory_layer_safety_2 = safety_set.get_current_memory_layer(DJ_safe);
					memory_layer_offset_check2 = seq_offsets.get_current_memory_layer(D_gene_seq,Three_prime);
					//cout<<"j_del_2: "<<memory_layer_safety_2<<endl;
				}
				break;
			default:
				break;
			}


	 //Get V 3' deletion
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime)) != 0){
		const Rec_Event* del_v_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,V_gene,Three_prime));
		if(processed_events.count(del_v_p->get_name())!=0){
			v_3_min_del=0;
			v_3_max_del=0;
		}
		else{
			v_3_min_del =  del_v_p->get_len_max();
			v_3_max_del =  del_v_p->get_len_min();
		}
	}
	else{
		v_3_min_del=0;
		v_3_max_del=0;
	}

	//Get D 5' deletion range
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime)) != 0){
		const Rec_Event* del_d_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Five_prime));
		if(processed_events.count(del_d_p->get_name())!=0){
			d_5_min_del=0;
			d_5_max_del=0;
		}
		else{
			d_5_min_del =  del_d_p->get_len_max();
			d_5_max_del =  del_d_p->get_len_min();
		}
	}
	else{
		d_5_min_del=0;
		d_5_max_del=0;
	}

	//Get D 3' deletion
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime)) != 0){
		const Rec_Event* del_d_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,D_gene,Three_prime));
		if(processed_events.count(del_d_p->get_name())!=0){
			d_3_min_del=0;
			d_3_max_del=0;
		}
		else{
			d_3_min_del =  del_d_p->get_len_max();
			d_3_max_del =  del_d_p->get_len_min();
		}
	}
	else{
		d_3_min_del=0;
		d_3_max_del=0;
	}

	//Get J 5' deletion range
	if(events_map.count(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime)) != 0){
		const Rec_Event* del_j_p = events_map.at(tuple<Event_type,Gene_class,Seq_side>(Deletion_t,J_gene,Five_prime));
		if(processed_events.count(del_j_p->get_name())!=0){
			j_5_min_del=0;
			j_5_max_del=0;
		}
		else{
			j_5_min_del= del_j_p->get_len_max();
			j_5_max_del= del_j_p->get_len_min();
		}
	}
	else{
		j_5_min_del=0;
		j_5_max_del=0;
	}
	this->Rec_Event::initialize_event(processed_events,events_map,offset_map,constructed_sequences,safety_set,error_rate_p,mismatches_list,seq_offsets,index_map);

 }

 void Deletion::add_to_marginals(long double scenario_proba , Marginal_array_p updated_marginals) const{
 	updated_marginals[this->new_index]+=scenario_proba;
 }


 string& make_transversions(string& init_sequence , bool is_int_seq){
	if(is_int_seq){
		for(string::iterator iter = init_sequence.begin() ; iter != init_sequence.end() ; ++iter){
			 if((*iter)=='0'){
				 (*iter) = '3';
			 }
			 else if ((*iter)=='1'){
				 (*iter) = '2';
			 }
			 else if ((*iter)=='2'){
				 (*iter) = '1';
			 }
			 else if ((*iter)=='3'){
				 (*iter)='0';
			 }
			 else if ((*iter)=='4'){
				 (*iter)='5';
			 }
			 else if ((*iter)=='5'){
				 (*iter)='4';
			 }
			 else if ((*iter)=='8'){
				 //Nothing to do
			 }
			 else if ((*iter)=='9'){
				 //Nothing to do
			 }
			 else if ((*iter)=='14'){
				 //Nothing to do
			 }
			 else{
				 throw runtime_error("Unknown int nucleotide " + to_string((*iter)) + " in seq " + init_sequence+" in make_transversions()");
			 }
		 }
	}
	else{
	 for(string::iterator iter = init_sequence.begin() ; iter != init_sequence.end() ; ++iter){
			 if((*iter)=='A'){
				 (*iter) = 'T';
			 }
			 else if ((*iter)=='C'){
				 (*iter) = 'G';
			 }
			 else if ((*iter)=='G'){
				 (*iter) = 'C';
			 }
			 else if ((*iter)=='T'){
				 (*iter)='A';
			 }
			 else{
				 throw runtime_error("Unknown int nucleotide " + to_string((*iter)) + " in seq " + init_sequence+" in make_transversions()");
			 }
		 }
	 }
	 return init_sequence;
 }


 bool del_numb_compare(const Event_realization& real1 , const Event_realization& real2) {
	 return real1.value_int > real2.value_int;
 }
