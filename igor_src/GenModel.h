/*
 * GenModel.h
 *
 *  Created on: 3 nov. 2014
 *      Author: marcou
 *      This class designs a generative model and supply all the methods to run a maximum likelihood estimate of the generative model
 */

#ifndef GENMODEL_H_
#define GENMODEL_H_

#include "Model_Parms.h"
#include "Rec_Event.h"
#include "Counter.h"
#include "Model_marginals.h"
#include "Errorrate.h"
#include "Utils.h"
#include <list>
#include <map>
#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <stdexcept>
#include <stack>
#include <memory>

//Make typedef for the function pointers
typedef void (*gen_seq_trans)(size_t , std::pair<std::string , std::queue<std::queue<int>>>,std::shared_ptr<void>);

/**
 * Hardcode a data structure for the function extracting CDR3s in generated sequences
 */
struct gen_CDR3_data{
	std::map<int,std::tuple<std::string,size_t,size_t,std::string>> v_anchors;
	size_t v_event_queue_position;
	std::map<int,std::tuple<std::string,size_t,size_t,std::string>> j_anchors;
	size_t j_event_queue_position;
	std::shared_ptr<std::ostream> output_stream;
	//Some config booleans
	bool output_nt_CDR3 = true;
	bool output_anchors_found = true;
	bool output_inframe  = true;
	//FIXME
	//For now do not output aa CDR3 stats
	bool output_aa_CDR3 = false;
	bool output_productive = false;


	gen_CDR3_data(const std::unordered_map<std::string,size_t>& v_anchors_indices , const std::unordered_map < std::string, Event_realization >& v_reals, size_t v_event_pos,
			const std::unordered_map<std::string,size_t>& j_anchors_indices , const std::unordered_map < std::string, Event_realization >& j_reals, size_t j_event_pos,
			std::shared_ptr<std::ostream> output_stream_ptr = std::shared_ptr<std::ostream>(&std::cout,null_delete<std::ostream>())): v_event_queue_position(v_event_pos) , j_event_queue_position(j_event_pos) , output_stream(output_stream_ptr){


		//First get all V anchors
		this->v_anchors.clear();
		for(const std::pair<std::string,Event_realization> v_real : v_reals){
			size_t v_anchor_index;
			if(v_anchors_indices.count(v_real.second.name)>0){
				v_anchor_index = v_anchors_indices.at(v_real.second.name);
				v_anchors.emplace(v_real.second.index,std::make_tuple(v_real.second.name,v_anchor_index,v_real.second.value_str.size(),v_real.second.value_str.substr(v_anchor_index,3)));
			}
			else{
				v_anchor_index = 0;
				v_anchors.emplace(v_real.second.index,std::make_tuple(v_real.second.name,v_anchor_index,v_real.second.value_str.size(),""));
			}
			/*try{
				v_anchor_index = v_anchors_indices.at(v_real.name);
			}
			catch (std::exception& e) {
				std::cerr<<"Could not find "<<v_real.name<<" in the V genes anchors map"<<std::endl;
				throw e;
			}*/
			//v_anchors.emplace(v_real.second.index,std::make_tuple(v_real.second.name,v_anchor_index,v_real.second.value_str.size(),v_real.second.value_str.substr(v_anchor_index,3)));
		}

		//Now get all J anchors
		this->j_anchors.clear();
		for(const std::pair<std::string,Event_realization> j_real : j_reals){
			size_t j_anchor_index;
			/*try{
				j_anchor_index = j_anchors_indices.at(j_real.name);
			}
			catch (std::exception& e) {
				std::cerr<<"Could not find "<<j_real.name<<" in the J genes anchors map"<<std::endl;
				throw e;
			}*/
			if(j_anchors_indices.count(j_real.second.name)>0){
				j_anchor_index = j_anchors_indices.at(j_real.second.name);
				j_anchors.emplace(j_real.second.index,std::make_tuple(j_real.second.name,j_anchor_index,j_real.second.value_str.size(),j_real.second.value_str.substr(j_anchor_index,3)));
			}
			else{
				j_anchor_index = std::string::npos;
				j_anchors.emplace(j_real.second.index,std::make_tuple(j_real.second.name,j_anchor_index,j_real.second.value_str.size(),""));
			}
			//j_anchors.emplace(j_real.second.index,std::make_tuple(j_real.second.name,j_anchor_index,j_real.second.value_str.size(),j_real.second.value_str.substr(j_anchor_index,3)));
		}

		//Write output file header
		*output_stream.get()<<"seq_index";
		if(output_nt_CDR3){
			*output_stream.get()<<",nt_CDR3";
		}
		if(output_anchors_found){
			*output_stream.get()<<",anchors_found";
		}
		if(output_inframe){
			*output_stream.get()<<",is_inframe";
		}
		if(output_aa_CDR3){
			*output_stream.get()<<",aa_CDR3";
		}
		if(output_productive){
			*output_stream.get()<<",is_productive";
		}
		*output_stream.get()<<std::endl;
	}
};

class GenModel {
public:
	GenModel(const Model_Parms&);
	GenModel(const Model_Parms& , const Model_marginals&);
	GenModel(const Model_Parms& , const Model_marginals& , const std::map<size_t,std::shared_ptr<Counter>>&);
	//TODO: add all the necessary constructors: with just model_parms, with model_parms and marginals
	virtual ~GenModel();

	bool infer_model(const std::vector<std::tuple<int,std::string,std::unordered_map<Gene_class , std::vector<Alignment_data>>>>& sequences ,const  int iterations ,const std::string path, bool fast_iter , double likelihood_threshold=1e-25 , bool viterbi_like=false);
	bool infer_model(const std::vector<std::tuple<int,std::string,std::unordered_map<Gene_class , std::vector<Alignment_data>>>>& sequences ,const  int iterations ,const std::string path, bool fast_iter=true , double likelihood_threshold=1e-25 , double proba_threshold_factor=0.001 );
	bool infer_model(const std::vector<std::tuple<int,std::string,std::unordered_map<Gene_class , std::vector<Alignment_data>>>>& sequences ,const  int iterations ,const std::string path, bool fast_iter , double likelihood_threshold , bool viterbi_like , double proba_threshold_factor , double mean_number_seq_err_thresh = INFINITY);

	std::forward_list<std::pair<std::string , std::queue<std::queue<int>>>> generate_sequences (int,bool);
	void generate_sequences(int,bool,std::string,std::string,std::list<std::pair<gen_seq_trans,std::shared_ptr<void>>> = std::list<std::pair<gen_seq_trans,std::shared_ptr<void>>>(),bool output_only_func = false , int=-1);
	bool load_genmodel();
	bool write2txt ();
	bool readtxt ();
	void write_seq2txt(std::string,std::forward_list<std::string>);
	void write_seq_real2txt(std::string , std::string , std::forward_list<std::pair<std::string , std::queue<std::queue<int>>>>);

	//write alignments, load alignments

private:
	Model_Parms model_parms;
	Model_marginals model_marginals;
	std::map<size_t,std::shared_ptr<Counter>> counters_list;//Size_t is a unique identifier for the Counter(useful for adding them up)
	std::pair<std::string , std::queue<std::queue<int>>> generate_unique_sequence(std::queue<std::shared_ptr<Rec_Event>> , std::unordered_map<Rec_Event_name,int> , const std::unordered_map<Rec_Event_name,std::vector<std::pair<std::shared_ptr<const Rec_Event>,int>>>& , std::default_random_engine& , bool =true);
	Model_marginals compute_marginals(std::list<std::string> sequences);
	Model_marginals compute_seq_marginals (std::string sequence);
	Model_marginals compute_seq_marginals (std::string sequence , std::list<std::list<std::string> > allowed_scenarios );

};

std::vector<std::tuple<int,std::string,std::unordered_map<Gene_class , std::vector<Alignment_data>>>> get_best_aligns (const std::vector<std::tuple<int,std::string,std::unordered_map<Gene_class , std::vector<Alignment_data>>>>&, Gene_class);



void output_CDR3_gen_data(size_t , std::pair<std::string , std::queue<std::queue<int>>> seq_and_real , std::shared_ptr<void> func_data);



#endif /* GENMODEL_H_ */
