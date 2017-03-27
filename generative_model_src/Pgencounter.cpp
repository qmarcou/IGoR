/*
 * Pgencounter.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: quentin
 */

#include "Pgencounter.h"

using namespace std;

Pgen_counter::Pgen_counter(): Pgen_counter("/tmp/" , false) {
}

Pgen_counter::Pgen_counter(std::string path): Pgen_counter(path , false) {
}

Pgen_counter::Pgen_counter(std::string path , bool output_Pgen_estimator_only , bool do_output_sequences): Counter(path) , output_sequences(do_output_sequences) , output_Pgen_estimator(output_Pgen_estimator_only) , sequence_Pgens_map(unordered_map<Int_Str,pair<double,long double>>()) , scenario_resulting_sequence(Int_Str()) , read_likelihood(0)  , v_gene(false) , d_gene(false) , j_gene(false) , vd_ins(false) , dj_ins(false) , vj_ins(false){
	if(output_Pgen_estimator and output_sequences){
		throw invalid_argument("Cannot set both \"output_Pgen_estimator_only\" and \"do_output_sequences\" to true. Pgen estimator is one line per read, otherwise every scenario sequence per read");
	}
	this->last_iter_only = true;
}


Pgen_counter::~Pgen_counter() {
	// TODO Auto-generated destructor stub
}

void Pgen_counter::initialize_counter(const Model_Parms& parms , const Model_marginals& marginals){
	if(not fstreams_created){
		output_pgen_file_ptr = shared_ptr<ofstream>(new ofstream);
		output_pgen_file_ptr->open(path_to_file + "Pgen_counts.csv");
		//Create the header
		if(output_Pgen_estimator){
			(*output_pgen_file_ptr.get())<<"seq_index;Pgen_estimate"<<endl;
		}
		else{
			if(output_sequences){
				(*output_pgen_file_ptr.get())<<"seq_index;scen_sequence;Pgen;P_joint_read_seq"<<endl;
			}
			else{
				(*output_pgen_file_ptr.get())<<"seq_index;Pgen;P_seq_given_read"<<endl;
			}
		}

		fstreams_created = true;
	}

	const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map = parms.get_events_map();
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
}

void Pgen_counter::count_scenario(long double scenario_seq_joint_proba , double scenario_probability , const string& original_sequence ,  Seq_type_str_p_map& constructed_sequences , const Seq_offsets_map& seq_offsets , const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map , Mismatch_vectors_map& mismatches_lists ){
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


	if (sequence_Pgens_map.count(scenario_resulting_sequence)>0){
		pair<double,long double>& Pgen_Pjoint_pair = sequence_Pgens_map[scenario_resulting_sequence];
		Pgen_Pjoint_pair.first+=scenario_probability;
		Pgen_Pjoint_pair.second+=scenario_seq_joint_proba;
	}
	else{
		pair<double,long double>& Pgen_Pjoint_pair = sequence_Pgens_map[scenario_resulting_sequence];
		//make proper initialization
		Pgen_Pjoint_pair.first=scenario_probability;
		Pgen_Pjoint_pair.second=scenario_seq_joint_proba;
	}

	read_likelihood+=scenario_seq_joint_proba;

}

void Pgen_counter::dump_sequence_data(int seq_index , int iteration_n ){

	double log_P_gen_estimate = 0;
	for(unordered_map<Int_Str,pair<double,long double>>::const_iterator iter = sequence_Pgens_map.begin() ; iter != sequence_Pgens_map.end() ; ++iter){
		if(output_Pgen_estimator){
			log_P_gen_estimate += (*iter).second.second/read_likelihood*log((*iter).second.first);
		}
		else{
			if(output_sequences){
				//(*output_pgen_file_ptr)<<seq_index<<";"<<(*iter).first<<";"<<(*iter).second.first<<";"<<(*iter).second.second/read_likelihood<<endl;
			}
			else{
				(*output_pgen_file_ptr.get())<<seq_index<<";"<<(*iter).second.first<<";"<<(*iter).second.second/read_likelihood<<endl;
			}
		}
	}
	if(output_Pgen_estimator){
		(*output_pgen_file_ptr.get())<<seq_index<<";"<<exp(log_P_gen_estimate)<<endl;
	}
	//Reset counters
	read_likelihood = 0.0;
	sequence_Pgens_map.clear();
}

void Pgen_counter::add_checked(shared_ptr<Counter> counter){
	return;
}

shared_ptr<Counter> Pgen_counter::copy() const{
	shared_ptr<Pgen_counter> counter_copy_ptr (new Pgen_counter());
	counter_copy_ptr->fstreams_created = this->fstreams_created;
	counter_copy_ptr->output_Pgen_estimator = this->output_Pgen_estimator;
	counter_copy_ptr->output_sequences = this->output_sequences;
	if(this->fstreams_created){
		counter_copy_ptr->output_pgen_file_ptr = this->output_pgen_file_ptr;
	}
	else{
		throw runtime_error("Counters should not be copied before stream initalization");
	}
	return counter_copy_ptr;
}
