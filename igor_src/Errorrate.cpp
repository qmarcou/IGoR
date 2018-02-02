/*
 * Errorrate.cpp
 *
 *  Created on: Jan 22, 2015
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
 */

#include "Errorrate.h"

using namespace std;

Error_rate::Error_rate():  debug_number_scenarios(0), model_log_likelihood(0) , number_seq(0) , seq_likelihood(0)  , seq_mean_error_number(0) , seq_probability(0) , updated(true) {
	// TODO Auto-generated constructor stub
	this->max_err = 0;
	this->max_noerr = 0;
}

Error_rate::~Error_rate() {
	// TODO Auto-generated destructor stub
}


void Error_rate::initialize(const unordered_map<tuple<Event_type,Gene_class,Seq_side>, shared_ptr<Rec_Event>>& events_map){
	//Do nothing
	//This method is called if no other method is supplied in the instantiated class
}

void Error_rate::norm_weights_by_seq_likelihood(Marginal_array_p& single_seq_marginal_array , const size_t marginal_array_size , const double seq_weight/*=1 by default*/){
	if(seq_likelihood!=0){
		for(size_t i = 0 ; i != marginal_array_size ; ++i){
			single_seq_marginal_array[i]/=this->seq_likelihood*seq_weight;
		}
	}
	else{
		//Everything should already be 0, just a debug check
		for(size_t i = 0 ; i != marginal_array_size ; ++i){
			single_seq_marginal_array[i]=0;
		}
	}
}

double Error_rate::get_seq_mean_error_number() const{
	if(seq_likelihood!=0){
		return seq_mean_error_number/seq_likelihood;
	}
	else{
		return 0;
	}
}



void add_to_err_rate(Error_rate* err_p1 , Error_rate* err_p2){
	if(err_p1->type() != err_p2->type()){
		throw invalid_argument("Cannot add error_rate of type " + err_p1->type() +" and " + err_p2->type() );
	}
	else{
		err_p1->add_checked(err_p2);
		return;
	}
}

