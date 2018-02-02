/*
 * Model_marginals.cpp
 *
 *  Created on: 3 nov. 2014
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

#include "Model_marginals.h"


using namespace std;

Model_marginals::Model_marginals() {
	marginal_array_smart_p = NULL;
	marginal_arr_size = -1;
}

Model_marginals::Model_marginals(const Model_Parms& model_parms) {
	marginal_arr_size = compute_size(model_parms);
	if(marginal_arr_size == 0){
		throw runtime_error("provided Model_parms imply empty marginals in Model_marginals::Model_marginals(const Model_Parms& model_parms)");
	}
	marginal_array_smart_p = Marginal_array_p(new long double[marginal_arr_size]);
	this->null_initialize();
}
/*
 * Provides a deep copy of the object
 */
Model_marginals::Model_marginals(const Model_marginals& other):marginal_array_smart_p(Marginal_array_p(new long double[other.marginal_arr_size])), marginal_arr_size(other.marginal_arr_size) {
	for(size_t i=0 ; i!=marginal_arr_size ; ++i){
		marginal_array_smart_p[i] = other.marginal_array_smart_p[i];
	}
}

Model_marginals::Model_marginals(size_t arr_size):marginal_arr_size(arr_size){
	marginal_array_smart_p = Marginal_array_p(new long double[marginal_arr_size]);
	this->null_initialize();
}

Model_marginals Model_marginals::empty_copy(){
	Model_marginals null_copy (this->marginal_arr_size);
	null_copy.debug_marg_name = "null tmp copy";
	return null_copy;
}

Model_marginals& Model_marginals::invert_edge(Rec_Event_name ev1_name , Rec_Event_name ev2_name , Model_Parms& model_parms){

	shared_ptr<Rec_Event> parent_ptr;
	shared_ptr<Rec_Event> child_ptr;

	/*
	 * Note here we touch the limit of the design of this code...
	 * There is no good way of checking that those marginals and those parms are related...
	 * TODO partially redesign?
	 */
	if(model_parms.has_edge(ev1_name,ev2_name)){
		parent_ptr = model_parms.get_event_pointer(ev1_name);
		child_ptr = model_parms.get_event_pointer(ev2_name);
	}
	else if(model_parms.has_edge(ev2_name,ev1_name)){
		parent_ptr = model_parms.get_event_pointer(ev2_name);
		child_ptr = model_parms.get_event_pointer(ev1_name);
	}
	else{
		throw runtime_error("Model_marginals::invert_edge() : the edge between \"" + ev1_name + "\" and \"" + ev2_name + "\" does not exist");
	}

	const unordered_map<Rec_Event_name,int> orig_marginals_index_map = this->get_index_map(model_parms);
	const unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> orig_marginals_inverse_offset_map = this->get_inverse_offset_map(model_parms);


	/*
	 * Make sure inverting the edge will not create a cycle by testing on a copy
	 */
		Model_Parms model_parms_test_copy( model_parms);
		try{
			model_parms_test_copy.invert_edge(ev1_name,ev2_name);
		}
		catch(exception& e){
			cerr<<"Exception caught trying to invert an edge on a test model parms in Model_marginals::invert_edge, the operation is most likely creating a cycle, throwing exception now..."<<endl;
			throw(e);
		}

	/*
	 * First recompute the marginals for the event losing a parent
	 * To do so we marginalize the joint probability of this event and the former parent over the former parent realizations
	 * This simply requires to marginalize the probability of the former parent over all events that are not parents of the former child
	 */
	list<shared_ptr<Rec_Event>> child_parents =  model_parms.get_parents(child_ptr);
	set<Rec_Event_name> child_kept_dependencies;
	for(const shared_ptr<Rec_Event> event_ptr : child_parents){
		child_kept_dependencies.emplace(event_ptr->get_name());
	}
	child_kept_dependencies.emplace(child_ptr->get_name()); //Add the child in order to get an array of the correct size directly

	pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> marginalized_parent_probabilities = compute_event_marginal_probability(parent_ptr->get_name(),child_kept_dependencies,model_parms);

	/*cout<<"Model_marginals::invert_edge()"<<endl;
	size_t tmp_size = 1;
	for(pair<Rec_Event_name,size_t> tmp : marginalized_parent_probabilities.first){
		cout<<tmp.first<<",";
		tmp_size*=tmp.second;
	}
	cout<<endl<<"Marginalized_parent_probabilities:"<<endl;
	for(size_t i=0 ; i!=tmp_size ; ++i){
		cout<<marginalized_parent_probabilities.second.get()[i]<<",";
	}
	cout<<endl;*/


	size_t joint_array_size = this->get_event_size(child_ptr,model_parms);
	size_t new_array_size = joint_array_size/parent_ptr->size();
	shared_ptr<long double> joint_child_array(new long double [joint_array_size]);
	shared_ptr<long double> new_child_array(new long double [new_array_size]);
	int child_index = orig_marginals_index_map.at(child_ptr->get_name());

	//Init the joint proba array
	for(size_t i=0 ; i!= joint_array_size ; ++i){
		joint_child_array.get()[i] = this->marginal_array_smart_p[child_index+i];
	}


	//Get the array ordering
	//First copy a list containing the inverse offsets and sort it
	list<pair<shared_ptr<const Rec_Event>,int>> sorted_inv_offset_list = orig_marginals_inverse_offset_map.at(child_ptr->get_name());
	sorted_inv_offset_list.sort(inverse_offset_comparator());

	//Create the ordering list
	list<pair<Rec_Event_name,size_t>> child_dependencies_order_list;
	child_dependencies_order_list.emplace_back(child_ptr->get_name() , child_ptr->size());
	//Now add dimensions in the correct order
	for(const pair<shared_ptr<const Rec_Event>,int>& inv_offset : sorted_inv_offset_list){
		child_dependencies_order_list.emplace_back(inv_offset.first->get_name(),inv_offset.first->size());
	}


	//Align the arrays
	align_marginal_array(child_dependencies_order_list,marginalized_parent_probabilities);

	//Multiply to get the joint
	for(size_t i=0 ; i!= joint_array_size ; ++i){
		joint_child_array.get()[i] *= marginalized_parent_probabilities.second.get()[i];
	}

	//Now move the former parent as the last dimension and marginalize over it
	pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> tmp_pair = make_pair(child_dependencies_order_list,joint_child_array);
	swap_events_order(parent_ptr->get_name() , child_dependencies_order_list.back().first , tmp_pair);

	for(size_t i=0 ; i!=new_array_size ; ++i){
		new_child_array.get()[i] = 0.0;
	}

	for(size_t j=0 ; j!=joint_array_size ; ++j){
		new_child_array.get()[j%new_array_size] += joint_child_array.get()[j];
	}
	/*cout<<"Marginalized joint child:"<<endl;
	for(size_t j=0 ; j!=new_array_size ; ++j){
		cout<<new_child_array.get()[j]<<",";
	}
	cout<<endl;*/

	/*
	 * Now we recompute the marginals of the event gaining a parent
	 * To do so we need to marginalize the former child distribution over all unshared parents except the former parent
	 * Multiply by the the former parent proba to get the joint
	 * Divide by the former child distribution marginalized over all unshared parents (including the former parent) for the bayesian inversion
	 */

	list<shared_ptr<Rec_Event>> parent_parents =  model_parms.get_parents(parent_ptr);
	set<Rec_Event_name> parent_kept_dependencies;
	for(const shared_ptr<Rec_Event> event_ptr : parent_parents){
		parent_kept_dependencies.emplace(event_ptr->get_name());
	}
	pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> marginalized_child_probabilities = compute_event_marginal_probability(child_ptr->get_name(),parent_kept_dependencies,model_parms);
	//By removing the dependence on the former parent we lose one dimension and we thus add it by hand
	//Compute the total size
	size_t new_child_array_size = parent_ptr->size();
	for(pair<Rec_Event_name,size_t> ev_name_size : marginalized_child_probabilities.first){
		new_child_array_size*= ev_name_size.second;
	}
	shared_ptr<long double> marginalized_child_proba_expanded_arr(new long double [new_child_array_size]);
	for(size_t i=0 ; i!=new_child_array_size ; ++i){
		marginalized_child_proba_expanded_arr.get()[i] = marginalized_child_probabilities.second.get()[i%(new_child_array_size/parent_ptr->size())];
	}
	marginalized_child_probabilities.first.emplace_back(parent_ptr->get_name() , parent_ptr->size());
	marginalized_child_probabilities.second = marginalized_child_proba_expanded_arr;

	/*cout<<"Marginalized child expanded: "<<endl;
	for(size_t i=0 ; i!=new_child_array_size ; ++i){
		cout<<marginalized_child_proba_expanded_arr.get()[i]<<",";
	}
	cout<<endl;
	cout<<"Sanity check on the marginalized child expanded pointer transfer:"<<endl;
	for(size_t i=0 ; i!=new_child_array_size ; ++i){
		cout<<marginalized_child_probabilities.second.get()[i]<<",";
	}
	cout<<endl;*/

	//Now compute the marginal distribution keeping the parent dependence
	parent_kept_dependencies.emplace(parent_ptr->get_name());
	pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> marginalized_child_g_parent_probabilities = compute_event_marginal_probability(child_ptr->get_name(),parent_kept_dependencies,model_parms);

	/*cout<<"Marginalized child given parent: "<<endl;
	size_t tmp_sizzze = 1;
	for(auto zz : marginalized_child_g_parent_probabilities.first){
		tmp_sizzze*=zz.second;
	}
	for(size_t i=0 ; i!=tmp_sizzze ; ++i){
		cout<<marginalized_child_g_parent_probabilities.second.get()[i]<<",";
	}
	cout<<endl;*/

	//Now extract the former marginal values and expand the dimension with the former child
	shared_ptr<long double> new_parent_array (new long double [new_child_array_size]);
	int parent_index = orig_marginals_index_map.at(parent_ptr->get_name());
	//cout<<"New parent array:"<<endl;
	for(size_t i=0 ; i!=new_child_array_size ; ++i){
		new_parent_array.get()[i] = this->marginal_array_smart_p.get()[parent_index + i%(new_child_array_size/child_ptr->size())];
		//cout<<new_parent_array.get()[i]<<",";
	}
	//cout<<endl;

	//Get the corresponding order list
	//First copy a list containing the inverse offsets and sort it
	sorted_inv_offset_list = orig_marginals_inverse_offset_map.at(parent_ptr->get_name());
	sorted_inv_offset_list.sort(inverse_offset_comparator());

	//Create the ordering list
	list<pair<Rec_Event_name,size_t>> parents_dependencies_order_list;
	parents_dependencies_order_list.emplace_back(parent_ptr->get_name() , parent_ptr->size());
	//Now add dimensions in the correct order
	for(const pair<shared_ptr<const Rec_Event>,int>& inv_offset : sorted_inv_offset_list){
		parents_dependencies_order_list.emplace_back(inv_offset.first->get_name(),inv_offset.first->size());
	}
	//And add the extra dimension
	parents_dependencies_order_list.emplace_back(child_ptr->get_name() , child_ptr->size());


	//Order everything (according to the new ordering in model parms)
	model_parms.invert_edge(ev1_name,ev2_name);
	const unordered_map<Rec_Event_name,int> new_marginals_index_map = this->get_index_map(model_parms);
	const unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> new_marginals_inverse_offset_map = this->get_inverse_offset_map(model_parms);

	//Create the ordering list to do so
	list<pair<Rec_Event_name,size_t>> final_new_child_ordering_list;
	final_new_child_ordering_list.emplace_back(parent_ptr->get_name() , parent_ptr->size());
	//Now add dimensions in the correct order
	for(const pair<shared_ptr<const Rec_Event>,int>& inv_offset : sorted_inv_offset_list){
		final_new_child_ordering_list.emplace_back(inv_offset.first->get_name(),inv_offset.first->size());
	}
	final_new_child_ordering_list.emplace_back(child_ptr->get_name() , child_ptr->size());

	//Align everything
	align_marginal_array(final_new_child_ordering_list,marginalized_child_g_parent_probabilities);
	align_marginal_array(final_new_child_ordering_list,marginalized_child_probabilities);
	pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> tmp_pair_parent = make_pair(parents_dependencies_order_list,new_parent_array);
	align_marginal_array(final_new_child_ordering_list,tmp_pair_parent);

/*	cout<<"Sanity test:"<<endl;
	size_t bobby = 1;
	for(auto zzz: marginalized_child_g_parent_probabilities.first){
		cout<<zzz.first<<",";
		bobby*=zzz.second;
	}
	cout<<endl;
	for(size_t i=0 ; i!= bobby ; ++i){
		cout<<marginalized_child_g_parent_probabilities.second.get()[i]<<",";
	}
	cout<<endl;
	bobby = 1;
	for(auto zzz: marginalized_child_probabilities.first){
		cout<<zzz.first<<",";
		bobby*=zzz.second;
	}
	cout<<endl;
	for(size_t i=0 ; i!= bobby ; ++i){
		cout<<marginalized_child_probabilities.second.get()[i]<<",";
	}
	cout<<endl;
	bobby = 1;
	for(auto zzz: tmp_pair_parent.first){
		cout<<zzz.first<<",";
		bobby*=zzz.second;
	}
	cout<<endl;
	for(size_t i=0 ; i!= bobby ; ++i){
		cout<<tmp_pair_parent.second.get()[i]<<",";
	}
	cout<<endl;
*/
	//Finally compute the new marginal values
	//cout<<"New parent array final:"<<endl;
	for(size_t i=0 ; i!=new_child_array_size ; ++i){
		if(marginalized_child_probabilities.second.get()[i]!=0){
			new_parent_array.get()[i] *= marginalized_child_g_parent_probabilities.second.get()[i]/marginalized_child_probabilities.second.get()[i];
		}
		else{
			new_parent_array.get()[i] = 0.0;
		}
		//cout<<new_parent_array.get()[i]<<",";
	}
	//cout<<endl;


	//Create a full new marginal array and copy the values in the right place
	unique_ptr<long double []> new_marginal_array (new long double[ this->compute_size(model_parms)]);
	for(const shared_ptr<Rec_Event> event_ptr : model_parms.get_event_list()){
		if( (event_ptr->get_name() != ev1_name) and (event_ptr->get_name() != ev2_name)){
			size_t tmp_event_size = this->get_event_size(event_ptr,model_parms);
			int former_event_index = orig_marginals_index_map.at(event_ptr->get_name());
			int new_event_index = new_marginals_index_map.at(event_ptr->get_name());
			for(size_t i=0 ; i!=tmp_event_size ; ++i){
				new_marginal_array[new_event_index+i] = this->marginal_array_smart_p[former_event_index+i];
			}
		}
	}
	//Now copy the new parent values
	int new_parent_index = new_marginals_index_map.at(child_ptr->get_name());
	for(size_t i=0 ; i!=new_array_size ; ++i){
		new_marginal_array[new_parent_index+i] = new_child_array.get()[i];
	}

	//And finally copy the new child values
	int new_child_index = new_marginals_index_map.at(parent_ptr->get_name());
	for(size_t i=0 ; i!=new_child_array_size ; ++i){
		new_marginal_array[new_child_index+i] = new_parent_array.get()[i];
	}

	//Set the array and update the array size
	this->marginal_arr_size = this->compute_size(model_parms);
	this->marginal_array_smart_p.swap(new_marginal_array);

	return *this;
}

size_t Model_marginals::compute_size(const Model_Parms& model_parms){
	list<shared_ptr<Rec_Event>> events = model_parms.get_event_list();
		unordered_map<Rec_Event_name,Adjacency_list> edges = model_parms.get_edges();

		size_t array_size = 0;
		for(list<shared_ptr<Rec_Event>>::const_iterator iter = events.begin() ; iter!= events.end() ; ++iter){

			int event_size = this->get_event_size((*iter),model_parms);

			array_size+=event_size;
		}
		return array_size;
}


/*
 *
 */
size_t Model_marginals::get_event_size(shared_ptr<const Rec_Event> event_p , const Model_Parms& model_parms ) const{
	size_t event_size = event_p->size();
	list<shared_ptr<Rec_Event>> parents_list = model_parms.get_edges().at( event_p->get_name() ).parents;

	for(list<shared_ptr<Rec_Event>>::const_iterator jiter = parents_list.begin() ; jiter!=parents_list.end();++jiter){
		event_size*=(*(*jiter)).size();
	}

	return event_size;
}

/*
 * Compute the event marginal probability distribution (free of dependencies)
 * \bug /!\ This function assumes the marginals are normalized /!\
 * FIXME make sure the marginals are normalized? otherwise defined up to a multiplicative constant??
 *
 * Compute the event marginal probability by recursion
 *
 * should return a list of offsets (maybe it would be better to return the event sizes) and corresponding events
 */
pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> Model_marginals::compute_event_marginal_probability(Rec_Event_name event_name , const Model_Parms& model_parms ) const{
	return this->compute_event_marginal_probability(event_name , set<Rec_Event_name>() ,model_parms);
}

pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> Model_marginals::compute_event_marginal_probability(Rec_Event_name event_name , const set<Rec_Event_name>& kept_dependencies_list , const Model_Parms& model_parms ) const{
	const unordered_map<Rec_Event_name,int> index_map = this->get_index_map(model_parms);
	const unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> offset_map = this->get_offsets_map(model_parms);
	const unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> inverse_offset_map = this->get_inverse_offset_map(model_parms);
	return this->compute_event_marginal_probability(event_name,kept_dependencies_list,model_parms,index_map,offset_map,inverse_offset_map);
}

/**
 * This piece of code is quite dirty, i strongly apologize for it
 */
pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> Model_marginals::compute_event_marginal_probability(
		Rec_Event_name event_name ,
		const set<Rec_Event_name>& kept_dependencies_list ,
		const Model_Parms& model_parms ,const unordered_map<Rec_Event_name,int>&  index_map ,
		const unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>>&  offset_map ,
		const unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>>& inverse_offset_map) const{


	//First get the event pointer, size and index
	shared_ptr<Rec_Event> event_ptr = model_parms.get_event_pointer(event_name);
	size_t event_size = event_ptr->size();
	size_t event_index = index_map.at(event_name);


	//Compute the total new array size and instantiate the corresponding array
	size_t new_array_size = event_size;
	for(Rec_Event_name ev_name : kept_dependencies_list){
		if(ev_name != event_name){
			new_array_size*= model_parms.get_event_pointer(ev_name)->size();
		}
	}
	shared_ptr<long double> marginal_proba_ptr (new long double [new_array_size]);

	//Initialize a dependencies list to hold array dimensions names
	list<pair<Rec_Event_name,size_t>> dependencies_order_list;

	//Get the list of the event's parents
	const list<shared_ptr<Rec_Event>> parents_list = model_parms.get_parents(event_name);

/*	cout<<"Model_marginals::compute_event_marginal_probability()"<<endl;
	cout<<"Marginalized event: "<<event_name<<endl<<"Kept dependencies: ";
	for(Rec_Event_name ev_name : kept_dependencies_list){
		cout<<ev_name<<",";
	}
	cout<<endl;*/

	//Now compute the marginal probabilities of the event realizations
	if(parents_list.empty()){
		/*
		 * If the event has no parents then the probabilities contained on the array are already the marginal probabilities
		 * This condition must be met recursively by reaching a root of the graph (this is ensured by the acyclicity of the graph)
		 * However the array still needs to be copied several times to match the kept dependencies format
		 */

		//Simply copy the values in the array (the number of times required by the kept dependencies)
		for(size_t i = 0 ; i!=new_array_size ; ++i){
			marginal_proba_ptr.get()[i] = this->marginal_array_smart_p[event_index+(i%event_size)];
		}
		//Emplace Event order
		dependencies_order_list.emplace_back(event_name,event_size);
		for(Rec_Event_name ev_name : kept_dependencies_list){
			if(ev_name != event_name){
				dependencies_order_list.emplace_back(ev_name , model_parms.get_event_pointer(ev_name)->size());
			}
		}
	}
	else{
		//Compute event marginal array size(including kept dependencies) and create an array on which we compute the joint probabilities
		size_t marginal_event_size = event_size;

		//Make a copy of the kept dependencies and remove the parents and the event itself from it in order to avoid duplicate dimensions
		set<Rec_Event_name> tmp_kept_dependencies_utility = kept_dependencies_list;

		/*
		 * Here we get the list of dimensions in the correct order with respect to the existing marginals
		 * We use the inverse offset map for it and the fact that the event itself is always the first dimension
		 */
		//Add the event itself as a first array dimension (since it always is the first dimension according to the offset map construction)
		dependencies_order_list.emplace_back(event_name,event_size);

		//Remove potential appearance of the event in the kept dependencies
		if(tmp_kept_dependencies_utility.count(event_name)>0){
			tmp_kept_dependencies_utility.erase(event_name);
		}

		//Create a list containing all parents if they are not already in the kept dependencies (avoid duplicates)
		//Use the inverse offset map to know the correct dimension ordering
		//Disclaimer: this is going to be quite ugly and should be rewritten //TODO

		//First copy a list containing the inverse offsets and sort it
		list<pair<shared_ptr<const Rec_Event>,int>> sorted_inv_offset_list = inverse_offset_map.at(event_name);
		sorted_inv_offset_list.sort(inverse_offset_comparator());

		/*cout<<"Inverse offset list: ";
		for(const pair<shared_ptr<const Rec_Event>,int>& inv_offset : sorted_inv_offset_list){
			cout<<"("<<inv_offset.first->get_name()<<","<<inv_offset.second<<");";
		}
		cout<<endl;*/


		//Now add dimensions in the correct order
		for(const pair<shared_ptr<const Rec_Event>,int>& inv_offset : sorted_inv_offset_list){

			dependencies_order_list.emplace_back(inv_offset.first->get_name(),inv_offset.first->size());
			marginal_event_size*=inv_offset.first->size();

			if(tmp_kept_dependencies_utility.count(inv_offset.first->get_name())>0){
				tmp_kept_dependencies_utility.erase(inv_offset.first->get_name());
			}
		}


		//Now append kept dependencies that are not parents or the event itself as extra dimensions
		for(Rec_Event_name kept_dep_name : tmp_kept_dependencies_utility){
			dependencies_order_list.emplace_back(kept_dep_name,model_parms.get_event_pointer(kept_dep_name)->size());
			marginal_event_size*=model_parms.get_event_pointer(kept_dep_name)->size();
		}


		// Instantiate a large array to record the full joint proba with kept dependencies
		shared_ptr<long double> joint_proba_array (new long double [marginal_event_size]);

		//Compute the joint probabilities
		//First copy the conditionals from the marginals (again this assumes we have the dimensions in the right order, extra kept dependencies are just copies)
		size_t event_original_marginal_array_size = this->get_event_size(model_parms.get_event_pointer(event_name),model_parms);
		for(size_t i=0 ; i!=marginal_event_size ; ++i){
			joint_proba_array.get()[i] = this->marginal_array_smart_p[event_index+i%event_original_marginal_array_size]; //THIS IS ASSUMING AN ORDER FOR THE MARGINALS, NEED TO MAKE SURE IT IS CORRECT
			//cout<<joint_proba_array.get()[i]<<",";
		}
		//cout<<endl;

		//For each parent multiply by the marginal probability to obtain the joint
		/*
		 * Here we call recursively the compute_event_marginal_probability function on every parent
		 * We append the full list of parents (necessary if some joint ancestors exist, or if one parent is ancestor of the other) to the ket dependencies
		 * We also append the event itself: although it cannot be an ancestor of the parents (acyclic graph) it will force to output an array of the correct size
		 * Note that this might result in very big arrays if the graph is big, and might turn out to be inefficient //TODO compute overlap with ancestors and reexpand after marginalization
		 */
		set<Rec_Event_name> new_kept_dependencies_list  = kept_dependencies_list;
		if(new_kept_dependencies_list.count(event_name)==0){
			new_kept_dependencies_list.emplace(event_name);
		}
		for(shared_ptr<Rec_Event> parent_event : parents_list){
			if(new_kept_dependencies_list.count(parent_event->get_name())==0){
				new_kept_dependencies_list.emplace(parent_event->get_name());
			}
		}

		for(shared_ptr<Rec_Event> parent_event : parents_list){

			//Compute the marginal parent proba
			pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> parent_marginal_proba = compute_event_marginal_probability(parent_event->get_name(),new_kept_dependencies_list,model_parms,index_map,offset_map,inverse_offset_map);

			/*cout<<"Marginalized parent order: ";
			for(auto zzzz : parent_marginal_proba.first){
				cout<<zzzz.first<<",";
			}
			cout<<endl;
			cout<<"Marginalized parent order array: "<<endl;
			for(size_t i=0 ; i!=marginal_event_size ; ++i){
				cout<<parent_marginal_proba.second.get()[i]<<",";
			}
			cout<<endl;*/

			//Now align the obtained marginals to the reference joint marginals
			align_marginal_array(dependencies_order_list,parent_marginal_proba);

			/*cout<<"Marginalized parent name: "<<parent_event->get_name()<<endl;
			cout<<"Marginalized parent reorder: ";
			for(auto zzzz : parent_marginal_proba.first){
				cout<<zzzz.first<<",";
			}
			cout<<endl;
			cout<<"Marginalized parent reorder array: "<<endl;
			double test_sum = 0.0;
			for(size_t i=0 ; i!=marginal_event_size ; ++i){
				cout<<parent_marginal_proba.second.get()[i]<<",";
				test_sum+=parent_marginal_proba.second.get()[i];
			}
			cout<<endl<<"Marginalized parent sum : ";
			cout<<test_sum<<endl;

			cout<<"original parent marginal: "<<endl;
			for(size_t i=0 ; i!=this->get_event_size(parent_event,model_parms) ; ++i){
				cout<<this->marginal_array_smart_p[index_map.at(parent_event->get_name())+i]<<",";
			}
			cout<<endl;*/

			//Once aligned simply multiply term by term to obtain the joint
			for(size_t i=0 ; i!=marginal_event_size ; ++i){
				joint_proba_array.get()[i] *= parent_marginal_proba.second.get()[i]; //THIS IS ASSUMING AN ORDER FOR THE MARGINALS? NEED TO MAKE SURE IT IS CORRECT
				//cout<<joint_proba_array.get()[i]<<",";
			}
			//cout<<endl;

		}

		/*double joint_sum = 0.0;
		cout<<"Joint proba array: "<<endl;
		for(size_t i=0 ; i!=marginal_event_size ; ++i){
			cout<<joint_proba_array.get()[i]<<',';
			joint_sum+=joint_proba_array.get()[i];
		}
		cout<<endl;
		cout<<"JOINT PROBA SUM: "<<joint_sum<<endl;*/


		//Now compute the final marginalized array
		//Initialize the array
		for(size_t i =0 ; i!= new_array_size ; ++i){
			marginal_proba_ptr.get()[i] = 0.0;
		}

		//Now re-organize the joint in order to have all necessary dimensions before marginalized ones
		list<pair<Rec_Event_name,size_t>> final_dimensions_order_list;
		final_dimensions_order_list.emplace_back(event_name,event_size);//Put the considered event in first position
		for(Rec_Event_name ev_name : kept_dependencies_list){
			if(ev_name != event_name){
				final_dimensions_order_list.emplace_back(ev_name , model_parms.get_event_pointer(ev_name)->size());
			}
		}
		pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>> tmp_pair = make_pair(dependencies_order_list,joint_proba_array);
		align_marginal_array(final_dimensions_order_list,tmp_pair);
		//Now reassign the correct pointer since i has been changed in the function and make_pair made a copy of teh smart pointer
		joint_proba_array = tmp_pair.second;
		dependencies_order_list = tmp_pair.first;

		/*cout<<"Joint proba array reorder: "<<endl;
		for(size_t i=0 ; i!=marginal_event_size ; ++i){
			cout<<joint_proba_array.get()[i]<<',';
		}
		cout<<endl;*/

		//Now from the joint compute the marginal joint probability (sum all extra dimensions)
		for(size_t i =0 ; i!= marginal_event_size ; ++i){
			marginal_proba_ptr.get()[i%new_array_size]+=joint_proba_array.get()[i];
		}

		//And finally get the conditional
		shared_ptr<long double> summed_joint_arr (new long double [new_array_size/event_size]);
		for(size_t i=0 ; i!= new_array_size/event_size ; ++i){
			summed_joint_arr.get()[i] = 0.0;
		}
		//compute the sum over dependencies
		for(size_t i =0 ; i!= new_array_size ; ++i){
			summed_joint_arr.get()[i/event_size]+=marginal_proba_ptr.get()[i];
		}
		//renormalize to get the conditional
		for(size_t i =0 ; i!= new_array_size ; ++i){
			if(summed_joint_arr.get()[i/event_size]!=0){
				marginal_proba_ptr.get()[i]/=summed_joint_arr.get()[i/event_size];
			}
			else{
				marginal_proba_ptr.get()[i] = 0.0;
			}
		}
		/*cout<<"Final dimension order list: ";
		for(auto zzzz : final_dimensions_order_list){
			cout<<zzzz.first<<",";
		}
		cout<<endl;
		cout<<"Dependencies order list: ";
		for(auto zzzz : dependencies_order_list){
			cout<<zzzz.first<<",";
		}
		cout<<endl;*/

		dependencies_order_list = final_dimensions_order_list;

	}

	return make_pair(dependencies_order_list,marginal_proba_ptr);
}

/*
 * Just a utility function to recursively
 */

Model_marginals::~Model_marginals() {
	//cout<<debug_marg_name<<endl;
	//delete [] this->marginal_array_p;
	//cout<<"test"<<endl;
}

Model_marginals& Model_marginals::operator=(const Model_marginals& other){
	//delete [] this->marginal_array_p;
	this->marginal_arr_size = other.marginal_arr_size;
	this->marginal_array_smart_p = Marginal_array_p(new long double [this->marginal_arr_size]);
	for(size_t i=0 ; i!=marginal_arr_size ; ++i){
		marginal_array_smart_p[i] = other.marginal_array_smart_p[i];
	}
	return *this;
}

Model_marginals& Model_marginals::operator +=(Model_marginals marginals){
	if(this->marginal_arr_size != marginals.marginal_arr_size){
		throw invalid_argument("Model_marginals must have the same size in : Model_marginals::operator+=");
	}
	else{
		for(size_t i = 0 ; i!= this->marginal_arr_size ; ++i){
			this->marginal_array_smart_p[i]+=marginals.marginal_array_smart_p[i];
		}
	}
	return *this;
}

Model_marginals& Model_marginals::operator -=(Model_marginals marginals){
	if(this->marginal_arr_size != marginals.marginal_arr_size){
		throw invalid_argument("Model_marginals must have the same size in : Model_marginals::operator+=");
	}
	else{
		for(size_t i = 0 ; i!= this->marginal_arr_size ; ++i){
			this->marginal_array_smart_p[i]-=marginals.marginal_array_smart_p[i];
		}
	}
	return *this;
}

Model_marginals Model_marginals::operator +(Model_marginals marginals){
	Model_marginals temp = *this;
	return temp+=marginals;
}

Model_marginals Model_marginals::operator -(Model_marginals marginals){
	Model_marginals temp = *this;
	return temp-=marginals;
}


unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_inverse_offset_map(const Model_Parms& model_parms) const{
	queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
	return get_inverse_offset_map(model_parms,model_queue);
}


unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_inverse_offset_map(const Model_Parms& model_parms, queue<shared_ptr<Rec_Event>> model_queue) const{
	//Stack to keep track of the events processed and their order
	stack<shared_ptr<Rec_Event>>* model_stack_p = new stack<shared_ptr<Rec_Event>>;
	stack<shared_ptr<Rec_Event>> model_stack = *model_stack_p;
		unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> invert_offset_map = unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> ();
		while(! model_queue.empty()){

			if(!model_stack.empty()){
				shared_ptr<Rec_Event> current_event_point = model_queue.front();
				//We look for all events upstream (parents)
				unordered_map<Rec_Event_name,shared_ptr<Rec_Event>>* related_events_map_p = new unordered_map<Rec_Event_name,shared_ptr<Rec_Event>>; //stores all the related events
				unordered_map<Rec_Event_name,shared_ptr<Rec_Event>> related_events_map = *related_events_map_p;

				if(!model_parms.get_parents( current_event_point ).empty()){
						list<shared_ptr<Rec_Event>> event_parents = model_parms.get_parents(current_event_point);

					for(list<shared_ptr<Rec_Event>>::const_iterator iter=event_parents.begin() ; iter!= event_parents.end() ; ++iter){
						related_events_map.insert(make_pair((*iter)->get_name(),*iter));
					}

				}
				//Copy model stack to have a modifiable copy
				//TODO model stack might be more complicated than a simple list
				stack<shared_ptr<Rec_Event>> model_stack_copy = model_stack;
				invert_offset_map.emplace(current_event_point->get_name(), list<pair< shared_ptr<const Rec_Event>,int>>() );

				int offset = (*current_event_point).size();

				while(!model_stack_copy.empty()){
					if(related_events_map.count(model_stack_copy.top()->get_name())>0){
						invert_offset_map.at(current_event_point->get_name()).push_back(make_pair(model_stack_copy.top(),offset));
						offset*= (*model_stack_copy.top()).size();
					}
					model_stack_copy.pop();
				}

				delete related_events_map_p;
			}

			model_stack.push(model_queue.front());
			model_queue.pop();
		}
		delete model_stack_p;
		return invert_offset_map;
}



/*
 * This method gives a map of offsets according to each event.
 * Every time a given event is chosen it will influence the position
 * where the marginals for event downstream have to be added.
 * E.g: the choice of a V (say V3) gene will influence the position on the array where p(ins|V) must be written
 * regardless of the number of insertions.
 */
unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_offsets_map(const Model_Parms& model_parms) const{
	queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
	return get_offsets_map(model_parms,model_queue);
}

unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_offsets_map(const Model_Parms& model_parms,queue<shared_ptr<Rec_Event>> model_queue) const{

	unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> invert_offset_map = get_inverse_offset_map(model_parms,model_queue);

	// Inversion of the map
	unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> offset_map;// = *(new unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>>);

	for(unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>>::const_iterator iter = invert_offset_map.begin() ; iter != invert_offset_map.end() ; ++iter){
		for(list<pair<shared_ptr<const Rec_Event>,int>>::const_iterator jiter = (*iter).second.begin() ; jiter != (*iter).second.end() ; ++jiter ){
			/*if(offset_map.count( (*jiter).first->get_name() ) == 0){
				offset_map.emplace( (*jiter).first->get_name() , *(new list<pair<shared_ptr<const Rec_Event>,int>>) );
			}*/

			offset_map[ (*jiter).first->get_name() ].push_back(make_pair( model_parms.get_event_pointer( (*iter).first ),(*jiter).second));
		}
	}
	return offset_map;

}



unordered_map<Rec_Event_name,int> Model_marginals::get_index_map(const Model_Parms& model_parms) const{
	queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
	return get_index_map( model_parms , model_queue);
}

unordered_map<Rec_Event_name,int> Model_marginals::get_index_map(const Model_Parms& model_parms , queue<shared_ptr<Rec_Event>> model_queue) const{


	int index = 0;
	stack<shared_ptr<Rec_Event>>* model_stack_p = new stack<shared_ptr<Rec_Event>>;
	stack<shared_ptr<Rec_Event>> model_stack = *model_stack_p;
	unordered_map<Rec_Event_name,int> index_map =  unordered_map<Rec_Event_name,int> ();

	while(! model_queue.empty()){

		//The index added in the map have been computed at the previous iteration (is 0 if first event in the queue)
		index_map.insert(make_pair(model_queue.front()->get_name(),index));
		shared_ptr<Rec_Event> current_event_point = model_queue.front();
		int event_size = current_event_point->size();

		/*if(!model_stack.empty()){

						//We look for all events upstream (parents)
						unordered_map<Rec_Event_name,shared_ptr<Rec_Event>>* related_events_map_p = new unordered_map<Rec_Event_name,shared_ptr<Rec_Event>>; //stores all the related events
						unordered_map<Rec_Event_name,shared_ptr<Rec_Event>> related_events_map = *related_events_map_p;

						if(!model_parms.get_parents( current_event_point ).empty()){
								list<shared_ptr<Rec_Event>> event_parents = model_parms.get_parents(current_event_point);

							for(list<shared_ptr<Rec_Event>>::const_iterator iter=event_parents.end() ; iter!= event_parents.end() ; iter++){
								related_events_map.insert(make_pair((*iter)->get_name(),*iter));
							}

						}
						//Copy model stack to have a modifiable copy
						//TODO model stack might be more complicated than a simple list
						stack<shared_ptr<Rec_Event>> model_stack_copy = model_stack;


				//Compute the event size: total size needed on the array to store informations on an event and the ones that it depends on
				while(!model_stack_copy.empty()){
					if(related_events_map.count(model_stack_copy.top()->get_name())!=0){
						event_size *= model_stack_copy.top()->size();
					}
					model_stack_copy.pop();
				}
				delete related_events_map_p;
		}*/
			list<shared_ptr<Rec_Event>> parents_list = model_parms.get_parents(current_event_point);
			if(!parents_list.empty()){
				for(list<shared_ptr<Rec_Event>>::const_iterator iter = parents_list.begin() ; iter != parents_list.end() ; ++iter ){
					event_size*=(*iter)->size();
				}
			}


			//This gives the pointer for the beginning of the next event(thus it will be inserted in the map at the next iteration)
			index+=event_size;
			model_stack.push(model_queue.front());
			model_queue.pop();

		}
	delete model_stack_p;

	return index_map;
}
/*
 * This method normalizes the marginal array so that each probability sums to 1.
 */
void Model_marginals::normalize(unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> inverse_offset_map , unordered_map<Rec_Event_name,int> index_map , queue<shared_ptr<Rec_Event>> model_queue){
	while (! model_queue.empty()){
		shared_ptr<Rec_Event> current_event_point = model_queue.front();
		list<pair< shared_ptr<const Rec_Event>,int>> related_events;
		if(inverse_offset_map.count(current_event_point->get_name()) == 0){
			related_events = list<pair<shared_ptr<const Rec_Event> , int >> (); //TODO change this to prevent memory leak
		}
		else{
			related_events = inverse_offset_map.at(current_event_point->get_name());
		}

		//TODO if list not empty
		iterate_normalize(current_event_point,related_events,index_map.at(current_event_point->get_name()) , 0 );
		//delete &related_events;

		model_queue.pop();
	}
}

void Model_marginals::iterate_normalize(shared_ptr<const Rec_Event> current_event_point, list<pair<shared_ptr<const Rec_Event>,int>>& related_events, int index , int current_offset){

	if(related_events.empty()) {
		/*double sum_marginals = 0;
		for(int iter =0 ; iter != current_event_point->size() ; iter++){
			sum_marginals+= this->marginal_array_p[index + current_offset + iter];
		}
		if(sum_marginals!=0){
			for(int iter =0 ; iter != current_event_point->size() ; iter++){
				this->marginal_array_p[index + current_offset + iter] /= sum_marginals;
			}
		}*/
		current_event_point->ind_normalize(this->marginal_array_smart_p,index + current_offset);


		return;
	}
	else{
		pair<shared_ptr<const Rec_Event>,int> processed_event = related_events.front();
		size_t processed_event_size = processed_event.first->size();
		int offset = processed_event.second;
		list<pair<shared_ptr<const Rec_Event>,int>> related_events_copy = related_events;
		related_events_copy.pop_front();
		for (size_t iter=0 ; iter!= processed_event_size ; ++iter){
			int new_offset = current_offset+iter*offset;
			iterate_normalize(current_event_point , related_events_copy , index ,new_offset); //TODO might be something fishy here
		}
	}
}

void Model_marginals::copy_fixed_events_marginals(const Model_marginals& source_marginals, const Model_Parms& parms,const unordered_map<Rec_Event_name,int>& index_map){
	list<shared_ptr<Rec_Event>> events = parms.get_event_list();
	for(list<shared_ptr<Rec_Event>>::const_iterator iter = events.begin() ; iter!=events.end() ; ++iter ){
		if((*iter)->is_fixed()){
			size_t event_size = this->get_event_size((*iter),parms);
			size_t first_index = index_map.at((*iter)->get_name());
			for(size_t i = first_index ; i!=first_index+event_size ; ++i){
				this->marginal_array_smart_p[i] = source_marginals.marginal_array_smart_p[i];
			}
		}
	}
}

/*
 * This method initializes the marginal array with uniform probability for each event.
 *
 */
void Model_marginals::uniform_initialize(const Model_Parms& parms){


	list<shared_ptr<Rec_Event>> events = parms.get_event_list();
	queue<shared_ptr<Rec_Event>> model_queue = parms.get_model_queue();
	unordered_map<Rec_Event_name,int> index_map = get_index_map(parms,model_queue);
	unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> inverse_offset_map = get_inverse_offset_map(parms,model_queue);


	for(size_t i=0 ; i!= compute_size(parms) ; ++i){
		marginal_array_smart_p[i]=1;
	}


	this->normalize(inverse_offset_map , index_map , model_queue);

}

void Model_marginals::null_initialize(){
	for(size_t i=0 ; i!= marginal_arr_size ; ++i){
			marginal_array_smart_p[i]=0;
		}
}

void Model_marginals::random_initialize(const Model_Parms& parms){
	//Create seed for random generator
	//create a seed from timer
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point time = myclock::now();
	myclock::duration dur = myclock::time_point::max() - time;

	unsigned time_seed = dur.count();
	//Instantiate random number generator
	default_random_engine generator =  default_random_engine(time_seed);
	uniform_real_distribution<double> distribution(0.0,1.0);
	for(size_t i = 0 ; i != this->marginal_arr_size ; ++i){
		marginal_array_smart_p[i] = distribution(generator);
	}

	list<shared_ptr<Rec_Event>> events = parms.get_event_list();
	queue<shared_ptr<Rec_Event>> model_queue = parms.get_model_queue();
	unordered_map<Rec_Event_name,int> index_map = get_index_map(parms,model_queue);
	unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> inverse_offset_map = get_inverse_offset_map(parms,model_queue);


	this->normalize(inverse_offset_map , index_map , model_queue);

}

void Model_marginals::flatten(shared_ptr<const Rec_Event> event,const Model_Parms& parms){
	queue<shared_ptr<Rec_Event>> model_queue = parms.get_model_queue();
	unordered_map<Rec_Event_name,int> index_map = get_index_map(parms,model_queue);
	unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> inverse_offset_map = get_inverse_offset_map(parms,model_queue);

	size_t event_size = this->get_event_size(event,parms);

	int base_index = index_map.at(event->get_name());

	//Set all values to 1 on the marginal array
	for(int i = base_index ; i!=base_index+event_size;++i){
		this->marginal_array_smart_p[i]=1;
	}
	this->normalize(inverse_offset_map,index_map,model_queue);

}

/**
 * Sets the realization probability to the given value
 * Note that the value will be set for all conditional dependences
 *
 * //TODO recode this in order to be able to set several realizations probas at the same time
 * //FIXME if the supplied new value is 1.0 there will be a zero division issue (this could be fixed by fixing all others to 0 instead of trying to set the supplied one to 1)
 */
void Model_marginals::set_realization_proba(string realization_name ,shared_ptr<const Rec_Event> event_ptr ,double new_value ,const Model_Parms& model_parms){
	if( (new_value<0) or (new_value>1) ){
		throw runtime_error("Invalid new probability value \"" + to_string(new_value) +"\" in Model_marginals::set_realization_proba()");
	}

	//First get the index_map and event marginal size
	const unordered_map<Rec_Event_name,int> index_map = this->get_index_map(model_parms);
	size_t marginal_event_size = this->get_event_size(event_ptr,model_parms);
	size_t event_size = event_ptr->size();

	//Get the realization index
	if(event_ptr->get_realizations_map().count(realization_name)<=0){
		throw runtime_error("Unknown realization \"" + realization_name
				+ "\" for event " + event_ptr->get_name() + " in Model_marginals::set_realization_proba");
	}
	const Event_realization event_real = event_ptr->get_realizations_map().at(realization_name);
	const size_t& real_index = event_real.index;
	const size_t& event_index = index_map.at(event_ptr->get_name());

	//Get the summed probabilities for all othe realizations for every conditioning
	double* summed_probas = new double[marginal_event_size/event_size];
	for(size_t i=0 ; i!=marginal_event_size ; ++i){
		if(i%event_size != real_index){
			summed_probas[i/event_size]+=this->marginal_array_smart_p[event_index+i];
		}
	}

	//Now compute the new value to set before normalization
	for(size_t i=0 ; i!=marginal_event_size/event_size ; ++i){
		summed_probas[i] *= new_value/(1-new_value);
	}

	//Finally set the value
	for(size_t i=0 ; i!=marginal_event_size/event_size ; ++i){
		this->marginal_array_smart_p[event_index+i*event_size+real_index] = summed_probas[i];
	}

	//Delete the array
	delete [] summed_probas;

	//Now renormalize the marginals
	auto inverse_offset_map = this->get_inverse_offset_map(model_parms);
	list<pair< shared_ptr<const Rec_Event>,int>> related_events;
	if(inverse_offset_map.count(event_ptr->get_name()) == 0){
		related_events = list<pair<shared_ptr<const Rec_Event> , int >> (); //TODO change this to prevent memory leak
	}
	else{
		related_events = inverse_offset_map.at(event_ptr->get_name());
	}

	this->iterate_normalize(event_ptr,related_events,event_index , 0 );
	return;
}

void Model_marginals::write2txt(string filename , const Model_Parms& model_parms){
	ofstream outfile(filename);
	queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
	unordered_map<Rec_Event_name,int> rank_map;
	list<shared_ptr<Rec_Event>> processed_events;
	unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> inv_offset_map = get_inverse_offset_map(model_parms , model_queue);
	unordered_map<Rec_Event_name,int> index_map = get_index_map(model_parms,model_queue);
	//Write into the file according to the model queue order
	while(!model_queue.empty()){
		shared_ptr<Rec_Event> current_event_p = model_queue.front();
		outfile<<"@"<<current_event_p->get_nickname()<<endl;
		outfile<<"$Dim[";
		//if(!processed_events.empty()){
			list<pair<shared_ptr<const Rec_Event>,int>> inv_offset_list = inv_offset_map[current_event_p->get_name()];
			if(!inv_offset_list.empty()){
				inv_offset_list.sort(offset_comp());
				for(list<pair<shared_ptr<const Rec_Event>,int>>::const_iterator jiter=inv_offset_list.begin() ; jiter != inv_offset_list.end() ; ++jiter){
					outfile<<(*jiter).first->size()<<",";
				}
			}
		//}
		outfile<<current_event_p->size()<<"]"<<endl;
		list<string> empty_str_list = list<string>();


		write2txt_iteration(inv_offset_list.begin() , inv_offset_list.end() , index_map[current_event_p->get_name()] , outfile , current_event_p , empty_str_list);
		model_queue.pop();


	}
}

void Model_marginals::write2txt_iteration(list<pair<shared_ptr<const Rec_Event>,int>>::const_iterator iter,const list<pair<shared_ptr<const Rec_Event>,int>>::const_iterator iter_end,int index,ofstream& outfile , shared_ptr<Rec_Event> current_event_p , list<string>& header){

	if(iter!=iter_end){
		for(int i = 0 ; i != (*iter).first->size() ; ++i){
			list<string> header_copy = header;
			header_copy.push_back(string("[")+(*iter).first->get_nickname()+string(",")+to_string(i)+string("]"));
			int new_index = index + i*(*iter).second;
			list<pair<shared_ptr<const Rec_Event>,int>>::const_iterator iter_copy = iter;
			++iter_copy;
			write2txt_iteration(iter_copy , iter_end , new_index , outfile , current_event_p , header_copy);
		}
	}
	else{
		outfile<<"#";
		if(!header.empty()){
			outfile<<header.front();
			list<string>::const_iterator base_jiter = header.begin();
			++base_jiter;

			for(list<string>::const_iterator jiter = base_jiter ; jiter != header.end() ; ++jiter){
				outfile<<","<<(*jiter);

			}

		}

		outfile<<endl;
		outfile<<"%"<<marginal_array_smart_p[index];
		for(int j=1 ; j < current_event_p->size() ; ++j){
			outfile<<","<<marginal_array_smart_p[index+j];
		}
		outfile<<endl;
	}
}

void Model_marginals::txt2marginals(string filename, const Model_Parms& parms){
	ifstream testfilestream(filename);
	if(!testfilestream){
		throw runtime_error("Unknown file: "+filename);
	}
	string line_str;
	//First count the marginals' size
	int size_counter = 0;
	while(getline(testfilestream,line_str)){
		if(line_str[0]=='%'){
			size_t semicolon_index =  line_str.find(",");
			++size_counter;
			while(semicolon_index!=string::npos){
				size_t next_comma_index = line_str.find(",", (semicolon_index+1) );
				semicolon_index = next_comma_index;
				++size_counter;
			}
		}
	}
	size_t current_marginals_size = this->compute_size(parms);
	if(size_counter!=current_marginals_size){
		throw runtime_error("Marginals contained in file \"" +filename+ "\" and supplied Model_Parms do not match in size. Make sure the Bayesian Network structure/event realizations  and the marginals are coherent.");
	}

	//Now read the actual marginals values
	ifstream infile(filename);
	int index = 0;
	while(getline(infile,line_str)){
		if(line_str[0]=='%'){
			size_t semicolon_index =  line_str.find(",");
			this->marginal_array_smart_p[index] = stod(line_str.substr(1,(semicolon_index)));
			++index;
			while(semicolon_index!=string::npos){
				size_t next_comma_index = line_str.find(",", (semicolon_index+1) );
				this->marginal_array_smart_p[index] = stod(line_str.substr( (semicolon_index+1) , (next_comma_index - semicolon_index -1) ));
				semicolon_index = next_comma_index;
				++index;
			}
		}
	}
	//this->marginal_arr_size = index;

	//Make sure marginals are normalized (deals with problem of float precision output from the text file)
	queue<shared_ptr<Rec_Event>> model_queue = parms.get_model_queue();
	unordered_map<Rec_Event_name,int> index_map = get_index_map(parms,model_queue);
	unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> inverse_offset_map = get_inverse_offset_map(parms,model_queue);

	this->normalize(inverse_offset_map , index_map , model_queue);
}

/**
 * A utility function to swap the order of the events on a marginal array (used to marginalize and invert edges)
 * This is used to further align the marginals and combine them
 */
void swap_events_order(const Rec_Event_name event_1 ,const Rec_Event_name event_2 , pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>>& swapped_marginals){
	//First get the positions of the events
	size_t event_1_position = 0;
	bool event_1_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_1_iterator;

	size_t event_2_position = 0;
	bool event_2_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_2_iterator;

//////////////////////////////////////
		/*cout<<"swap_events_order()"<<endl;
		cout<<event_1<<endl;
		cout<<event_2<<endl;
		size_t tmp_size = 1;
		for(pair<Rec_Event_name,size_t> tmp : swapped_marginals.first){
			cout<<tmp.first<<",";
			tmp_size*=tmp.second;
		}
		cout<<endl;
		for(size_t i=0 ; i!=tmp_size ; ++i){
			cout<<swapped_marginals.second.get()[i]<<",";
		}
		cout<<endl;*/
/////////////////////////////////////

	for(list<pair<Rec_Event_name,size_t>>::iterator iter = swapped_marginals.first.begin() ; iter!=swapped_marginals.first.end() ; ++iter){
		if(iter->first == event_1){
			event_1_found = true;
		}
		if(iter->first == event_2){
			event_2_found = true;
		}

		if(not event_1_found){
			++event_1_position;
		}
		if(not event_2_found){
			++event_2_position;
		}
	}

	//Make sure both events were found
	if(not (event_1_found and event_2_found)){
		throw runtime_error(event_1 + " and " + event_2 + " not found on the array, in swap_events_order()");
	}

	//In case event 1 and 2 were the same events, no change
	if(event_1_position == event_2_position){return;}

	//Then swap neighbors with the first event until it reaches the correct position
	int increment_factor = (event_1_position < event_2_position)*2 -1;
	size_t initial_event_1_position = event_1_position;
	while(event_1_position != event_2_position){
		for(list<pair<Rec_Event_name,size_t>>::iterator iter = swapped_marginals.first.begin() ; iter!=swapped_marginals.first.end() ; ++iter){
			if(iter->first == event_1){
				event_1_iterator = iter;
				break;
			}
		}
		list<pair<Rec_Event_name,size_t>>::iterator next_event_iter = event_2_iterator;
		if(increment_factor==-1){
			--next_event_iter;
		}
		else if(increment_factor==1){
			++next_event_iter;
		}
		else{
			throw runtime_error("issue in swap event order");
		}
		//swap_neighboring_events_order(event_1 , (event_1_iterator+increment_factor)->first , swapped_marginals);
		swap_neighboring_events_order(event_1 , next_event_iter->first , swapped_marginals);
		event_1_position+=increment_factor;
	}
	//Adjust event_2_position since during the last move it has been swapped with event_1
	event_2_position-=increment_factor;

	//Then swap down the second one
	increment_factor*=(-1); //Now change the swapping direction
	while(event_2_position != initial_event_1_position){
		for(list<pair<Rec_Event_name,size_t>>::iterator iter = swapped_marginals.first.begin() ; iter!=swapped_marginals.first.end() ; ++iter){
			if(iter->first == event_2){
				event_2_iterator = iter;
				break;
			}
		}
		list<pair<Rec_Event_name,size_t>>::iterator next_event_iter = event_1_iterator;
		if(increment_factor==-1){
			--next_event_iter;
		}
		else if(increment_factor==1){
			++next_event_iter;
		}
		else{
			throw runtime_error("issue in swap event order");
		}
		//swap_neighboring_events_order(event_1 , (event_2_iterator+increment_factor)->first , swapped_marginals);
		swap_neighboring_events_order(event_2 , next_event_iter->first , swapped_marginals);
		event_2_position+=increment_factor;
	}

	/*cout<<"swap_events_order() suite"<<endl;
	cout<<"Swapped version"<<endl;
	cout<<event_2<<endl;
	tmp_size = 1;
	for(pair<Rec_Event_name,size_t> tmp : swapped_marginals.first){
		cout<<tmp.first<<",";
		tmp_size*=tmp.second;
	}
	cout<<endl;
	for(size_t i=0 ; i!=tmp_size ; ++i){
		cout<<swapped_marginals.second.get()[i]<<",";
	}
	cout<<endl;*/
}
/**
 * \bug Will invalidate iterators to the swapped marginals
 */
void swap_neighboring_events_order(const Rec_Event_name event_1 ,const Rec_Event_name event_2 , pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>>& swapped_marginals){
	size_t event_1_offset;
	size_t event_1_new_offset;
	bool event_1_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_1_iterator;

	size_t event_2_offset;
	size_t event_2_new_offset;
	bool event_2_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_2_iterator;

	size_t current_offset = 1;

	/*cout<<"swap_neighboring_events_order : "<<endl;
	cout<<event_1<<endl;
	cout<<event_2<<endl;
	size_t tmp_size = 1;
	for(pair<Rec_Event_name,size_t> tmp : swapped_marginals.first){
		cout<<tmp.first<<",";
		tmp_size*=tmp.second;
	}
	cout<<endl;
	for(size_t i=0 ; i!=tmp_size ; ++i){
		cout<<swapped_marginals.second.get()[i]<<",";
	}
	cout<<endl;*/

	for(list<pair<Rec_Event_name,size_t>>::iterator iter = swapped_marginals.first.begin() ; iter!=swapped_marginals.first.end() ; ++iter){
		if( iter->first == event_1 and not event_2_found){
			event_1_offset = current_offset;
			event_1_iterator = iter;
			event_1_found = true;
			//Update the current offset accordingly
			//current_offset*=iter->second;

		}
		else if( iter->first == event_2 and not event_1_found){


			event_2_offset = current_offset;
			event_2_iterator = iter;
			event_2_found = true;
			//Update the current offset accordingly
			//current_offset*=iter->second;

		}

		current_offset*=iter->second;//Need to get the total array size, thus loop over all realizations
		//cout<<iter->first<<",";
	}
	//cout<<endl;

	if(event_1_found){
		//Event 1 is first on the array, event_2 should be right after

		event_2_iterator = event_1_iterator;
		++event_2_iterator;

		if(event_2_iterator->first != event_2){
			throw runtime_error(event_1 + " and " + event_2 + " are not neighbors, in swap_neighboring_events_order()");
		}
		event_1_new_offset =event_1_offset*event_2_iterator->second;
		event_2_offset = event_1_offset*event_1_iterator->second;
		event_2_new_offset = event_1_offset;
	}
	else if(event_2_found){
		//Event 2 is first on the array, event_1 should be right after

		event_1_iterator = event_2_iterator;
		++event_1_iterator;

		if(event_1_iterator->first != event_1){
			throw runtime_error(event_1 + " and " + event_2 + " are not neighbors, in swap_neighboring_events_order()");
		}
		event_2_new_offset =event_2_offset*event_1_iterator->second;
		event_1_offset = event_2_offset*event_2_iterator->second;
		event_1_new_offset = event_2_offset;
	}
	else{
		//Throw exception if none of the two events were found
		throw runtime_error(event_1 + " and " + event_2 + " not found on the array, in swap_neighboring_events_order()");
	}

	shared_ptr<long double> new_array_ptr(new long double [current_offset]); //Current offset is the total size of the array after looping over all events

	for(size_t i=0 ; i!=current_offset ; ++i){ //Current offset has been used to compute the total size of the array
		if(event_1_found){
			//event 1 used to have
		}
		size_t small_offset = min(event_1_new_offset,event_2_new_offset);
		size_t big_offset = max(event_1_new_offset*event_1_iterator->second,event_2_new_offset*event_2_iterator->second);
		/*
		 * Now invert the marginals
		 * (i/former_offset)%size corresponds to the realization index
		 * i%small offset are all lower dependencies that need to be copied in the same order (lower dimensions)
		 * big offset denotes the offset of stuff depending on the two events (higher dimensions)
		 *
		 */
		new_array_ptr.get()[ (i/big_offset)*big_offset +
					   ((i/event_1_offset)%event_1_iterator->second)*event_1_new_offset +
					   ((i/event_2_offset)%event_2_iterator->second)*event_2_new_offset +
					   i%small_offset ] = swapped_marginals.second.get()[i];
	}


	//The STL does not provide a function to swap two elements positions
	//Instead i'll insert the two elements in the swapped order after the 2 elements in the right order
	//I will then delete the two unswapped elements
/*	if(event_1_found){
	/*	//swapped_marginals.first.insert(event_2_iterator+1,event_1_iterator,event_2_iterator+1);
		//swapped_marginals.first.erase(event_1_iterator,event_2_iterator+1);
		list<pair<Rec_Event_name,size_t>>::iterator tmp_iter = event_2_iterator;
		++tmp_iter;
		swapped_marginals.first.insert(tmp_iter,event_1_iterator,tmp_iter);
		swapped_marginals.first.erase(event_1_iterator,tmp_iter);

	}
	else if(event_2_found){
		//swapped_marginals.first.insert(event_1_iterator+1,event_2_iterator,event_1_iterator+1);
		//swapped_marginals.first.erase(event_2_iterator,event_1_iterator+1);
		list<pair<Rec_Event_name,size_t>>::iterator tmp_iter = event_1_iterator;
		++tmp_iter;
		swapped_marginals.first.insert(tmp_iter,event_2_iterator,tmp_iter);
		swapped_marginals.first.erase(event_2_iterator,tmp_iter);
	}
	*/
	//Now swap the neighboring events in the list
	std::swap(*event_1_iterator,*event_2_iterator);

	//Reassign the swapped array
	swapped_marginals.second = new_array_ptr;

	/*cout<<"swap_neighboring_events_order swapped : "<<endl;
	tmp_size = 1;
	for(pair<Rec_Event_name,size_t> tmp : swapped_marginals.first){
		cout<<tmp.first<<",";
		tmp_size*=tmp.second;
	}
	cout<<endl;
	for(size_t i=0 ; i!=tmp_size ; ++i){
		cout<<swapped_marginals.second.get()[i]<<",";
	}
	cout<<endl;*/

}

/**
 * Utility to align marginals with events in the same order
 * Note that the implementation is probably not efficient but at least remains simple
 * Can align marginals with higher number of dimensions than the reference
 */
void align_marginal_array(const list<pair<Rec_Event_name,size_t>>& reference_marginals_order , pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<long double>>& aligned_marginals){

	////////////////////////////////////////////////////
	/*cout<<"align_marginal_array()"<<endl;
	cout<<"Reference marginal order: ";
	for(pair<Rec_Event_name,size_t> zzz : reference_marginals_order){
		cout<<zzz.first<<",";
	}
	cout<<endl<<"Aligned marginal order: ";
	for(pair<Rec_Event_name,size_t> zzz : aligned_marginals.first){
		cout<<zzz.first<<",";
	}
	cout<<endl;*/


	if(reference_marginals_order.size()>aligned_marginals.first.size()){
		throw runtime_error("Aligned marginals have less dimensions than the reference marginals in align_marginal_array()");
	}
	list<pair<Rec_Event_name,size_t>>::const_iterator reference_iterator = reference_marginals_order.cbegin();
	size_t counter = 0;
	while(reference_iterator != reference_marginals_order.cend()){
		//Swap the reference event directly at the right position
		list<pair<Rec_Event_name,size_t>>::iterator tmp_iter = aligned_marginals.first.begin();
		size_t i=0;
		while(i!=counter){
			++tmp_iter;
			++i;
		}

		if(reference_iterator->first != tmp_iter->first){
			swap_events_order(reference_iterator->first , tmp_iter->first , aligned_marginals);
		}

		++counter;
		++reference_iterator;
	}
}



