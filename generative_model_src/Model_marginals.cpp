/*
 * Model_marginals.cpp
 *
 *  Created on: 3 nov. 2014
 *      Author: marcou
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
	const unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> orig_marginals_offset_map = this->get_offsets_map(model_parms);





	//Compute the new marginals of the former child
	//For the former child simply need to marginalize over the former parent
	const size_t child_orig_index = orig_marginals_index_map.at(parent_ptr->get_name());
	const size_t child_size = child_ptr->size();
	const size_t parent_size = parent_ptr->size();
	const size_t child_former_marginal_size = this->get_event_size(child_ptr,model_parms);
	double child_parent_child_joint_proba_array [child_former_marginal_size];
	pair<size_t,shared_ptr<double>> parent_marginal_proba = compute_event_marginal_probability(parent_ptr->get_name(),model_parms);

	//Get the parent offset
	size_t parent_offset;
	const vector<pair<shared_ptr<const Rec_Event>,int>>& offset_vec = orig_marginals_offset_map.at(parent_ptr->get_name());
	for(pair<shared_ptr<const Rec_Event>,int> offset_pair : offset_vec){
		if(offset_pair.first == child_ptr){
			//Again only comparing pointer addresses instead of the actual events values, should work in principle
			parent_offset = offset_pair.second;
		}
	}
	//Compute the joint
	for(size_t i=0 ; i!= parent_size; ++i){
		for(size_t j=0 ; j!=parent_offset ; ++j){
			child_parent_child_joint_proba_array[i*parent_offset + j] = this->marginal_array_smart_p.get()[child_orig_index + i*parent_offset + j]*parent_marginal_proba.second.get()[i];
		}
	}

	//Compute the new marginals for the former parent
	//For the former parent need to compute the full joint
	size_t parent_orig_index = orig_marginals_index_map.at(parent_ptr->get_name());


	//Copy all unchanged values of the former marginals


	//Note: Would it matter to have unormalized marginals?


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
pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<double>> Model_marginals::compute_event_marginal_probability(Rec_Event_name event_name , const Model_Parms& model_parms ) const{
	return this->compute_event_marginal_probability(event_name , list<shared_ptr<Rec_Event>>() ,model_parms);
}

pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<double>> Model_marginals::compute_event_marginal_probability(Rec_Event_name event_name , const list<shared_ptr<Rec_Event>>& kept_dependencies_list , const Model_Parms& model_parms ) const{
	const unordered_map<Rec_Event_name,int> index_map = this->get_index_map(model_parms);
	const unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> offset_map = this->get_offsets_map(model_parms);
	return this->compute_event_marginal_probability(event_name,kept_dependencies_list,model_parms,index_map,offset_map);
}

pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<double>> Model_marginals::compute_event_marginal_probability(Rec_Event_name event_name , const list<shared_ptr<Rec_Event>>& kept_dependencies_list , const Model_Parms& model_parms ,const unordered_map<Rec_Event_name,int>&  index_map , const unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>>&  offset_map) const{
	shared_ptr<Rec_Event> event_ptr = model_parms.get_event_pointer(event_name);
	size_t event_size = event_ptr->size();
	//Compute the total new array size
	size_t new_array_size = event_size;
	for(shared_ptr<Rec_Event> event_ptr : kept_dependencies_list){
		new_array_size*=event_ptr->size();
	}
	shared_ptr<double> marginal_proba_ptr (new double [new_array_size]);
	list<pair<Rec_Event_name,size_t>> dependencies_order_list;
	size_t event_index = index_map.at(event_name);
	const list<shared_ptr<Rec_Event>> parents_list = model_parms.get_parents(event_name);

	//Now compute the marginal probabilities of the event realizations
	if(parents_list.empty()){
		/*
		 * If the event has no parents then the probabilities contained on the array are already the marginal probabilities
		 * This condition must be met recursively by reaching a root of the graph (this is ensured by the acyclicity of the graph)
		 * However the array still needs to be copied several times to match the kept dependencies format
		 */

		//Simply copy the values in the array (the number of times required by the kept dependencies)
		for(size_t i = 0 ; i!=event_size ; ++i){
			marginal_proba_ptr.get()[i] = this->marginal_array_smart_p.get()[event_index+(i%event_size)];
		}
		//Emplace Event order
		dependencies_order_list.emplace_back(event_name,event_size);
		for(shared_ptr<Rec_Event> event_ptr : kept_dependencies_list){
			dependencies_order_list.emplace_back(event_ptr->get_name() , event_ptr->size());
		}
	}
	else{
		//Compute event marginal array size and create an array on which we compute the joint probabilities
		size_t marginal_event_size = this->get_event_size(event_ptr,model_parms);

		//Create a utility list
		joint_proba_array
		//Preprocess parents and see the overlap with kept dependencies
		for(shared_ptr<Rec_Event> event_ptr : dependencies_order_list){

		}

		double joint_proba_array [marginal_arr_size];

		//Compute the joint probabilities
		//First copy the conditionals from the marginals
		for(size_t i=0 ; i!=marginal_event_size ; ++i){
			joint_proba_array[i] = this->marginal_array_smart_p[event_index+i];
		}
		//For each parent multiply by the marginal probability to obtain the joint
		for(shared_ptr<Rec_Event> parent_event : parents_list){
			pair<size_t,shared_ptr<double>> parent_marginal_proba = compute_event_marginal_probability(parent_event->get_name(),kept_dependencies_list,model_parms,index_map,offset_map);
			size_t parent_size = parent_event->size();
			size_t parent_offset;
			bool event_found = false; //Should always be found, just a sanity check
			const vector<pair<shared_ptr<const Rec_Event>,int>>& offsets_vector = offset_map.at(parent_event->get_name());
			for(pair<shared_ptr<const Rec_Event>,int> offset_pair : offsets_vector){
				if(offset_pair.first == event_ptr){
					parent_offset = offset_pair.second;
					event_found = true;
					/*
					 * Note: comparison is made on the pointer address and not the event itself
					 * This should not be a problem a priori, otherwise the following exception will probably be raised
					 */
				}
				if(event_found){
					for(size_t i=0 ; i!=parent_size ; ++i){
						for(size_t j=0 ; j!=parent_offset ; ++j){
							joint_proba_array[i*parent_offset + j]*=parent_marginal_proba.second.get()[i];
						}
					}
				}
				else{
					throw logic_error("If this exception is thrown there is a problem in the Model_marginals::compute_event_marginal_probability code");
				}
			}
		}

		//Marginalize
		//Initialize the array
		for(size_t i =0 ; i!= event_size ; ++i){
			marginal_proba_ptr.get()[i] = 0.0;
		}
		//Now from the joint compute the marginal probability
		for(size_t i =0 ; i!= marginal_event_size ; ++i){
			marginal_proba_ptr.get()[i%event_size]+=joint_proba_array[i];
		}
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
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("Unknown file: "+filename);
	}
	string line_str;
	int index = 0;
	while(getline(infile,line_str)){
		if(line_str[0]=='%'){
			size_t semicolon_index =  line_str.find(",");
			marginal_array_smart_p[index] = stod(line_str.substr(1,(semicolon_index)));
			++index;
			while(semicolon_index!=string::npos){
				size_t next_comma_index = line_str.find(",", (semicolon_index+1) );
				marginal_array_smart_p[index] = stod(line_str.substr( (semicolon_index+1) , (next_comma_index - semicolon_index -1) ));
				semicolon_index = next_comma_index;
				++index;
			}
		}
	}
	this->marginal_arr_size = index;

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
void swap_events_order(const Rec_Event_name& event_1 ,const Rec_Event_name& event_2 , pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<double>>& swapped_marginals){
	//First get the positions of the events
	size_t event_1_position = 0;
	bool event_1_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_1_iterator;

	size_t event_2_position = 0;
	bool event_2_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_2_iterator;


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
		if(not event_1_found){
			++event_2_position;
		}
	}

	if(not (event_1_found and event_2_found)){
		throw runtime_error(event_1 + " and " + event_2 + " not found on the array, in swap_events_order()");
	}
	if(event_1_position == event_2_position){return;} //In case event 1 and 2 were the same events, no change

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
		swap_neighboring_events_order(event_1 , (event_1_iterator+increment_factor)->first , swapped_marginals);
		event_1_position+=increment_factor;
	}

	//Then swap down the second one
	increment_factor*=(-1); //Now change the swapping direction
	while(event_2_position != initial_event_1_position){
		for(list<pair<Rec_Event_name,size_t>>::iterator iter = swapped_marginals.first.begin() ; iter!=swapped_marginals.first.end() ; ++iter){
			if(iter->first == event_2){
				event_2_iterator = iter;
				break;
			}
		}
		swap_neighboring_events_order(event_1 , (event_2_iterator+increment_factor)->first , swapped_marginals);
		event_2_position+=increment_factor;
	}
}
/**
 * \bug Will invalidate iterators to the swapped marginals
 */
void swap_neighboring_events_order(const Rec_Event_name& event_1 ,const Rec_Event_name& event_2 , pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<double>>& swapped_marginals){
	size_t event_1_offset;
	size_t event_1_new_offset;
	bool event_1_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_1_iterator;

	size_t event_2_offset;
	size_t event_2_new_offset;
	bool event_2_found = false;
	list<pair<Rec_Event_name,size_t>>::iterator event_2_iterator;

	size_t current_offset = 1;

	for(list<pair<Rec_Event_name,size_t>>::iterator iter = swapped_marginals.first.begin() ; iter!=swapped_marginals.first.end() ; ++iter){
		if( iter->first == event_1 and not event_2_found){
			event_1_offset = current_offset;
			event_1_iterator = iter;
			event_1_found = true;
			//Update the current offset accordingly
			current_offset*=iter->second;

		}
		else if( iter->first == event_2 and not event_1_found){


			event_2_offset = current_offset;
			event_2_iterator = iter;
			event_2_found = true;
			//Update the current offset accordingly
			current_offset*=iter->second;

		}

		current_offset*=iter->second;//Need to get the total array size, thus loop over all realizations

	}

	if(event_1_found){
		//Event 1 is first on the array

		event_2_iterator = event_1_iterator+1;
		if(event_2_iterator->first != event_2){
			throw runtime_error(event_1 + " and " + event_2 + " are not neighbors, in swap_neighboring_events_order()");
		}
		event_1_new_offset =event_1_offset*event_2_iterator->second;
		event_2_offset = event_1_offset*event_1_iterator->second;
		event_2_new_offset = event_1_offset;
	}
	else if(event_2_found){
		//Event 2 is first on the array

		event_1_iterator = event_2_iterator+1;
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

	shared_ptr<double> new_array_ptr(new double [current_offset]); //Current offset is the total size of the array after looping over all events

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
		new_array_ptr[ (i/big_offset)*big_offset +
					   ((i/event_1_offset)%event_1_iterator->second)*event_1_new_offset +
					   ((i/event_2_offset)%event_2_iterator->second)*event_2_new_offset +
					   i%small_offset ] = swapped_marginals.second[i];
	}

	//Now swap the neighboring events in the list
	//The STL does not provide a function to swap to elements positions
	//Instead i'll insert the two elements in the swapped order after the 2 elements in the right order
	//I will then delete the two unswapped elements
	if(event_1_found){
		swapped_marginals.first.insert(event_2_iterator+1,event_1_iterator,event_2_iterator+1);
		swapped_marginals.first.erase(event_1_iterator,event_2_iterator+1);
	}
	else if(event_2_found){
		swapped_marginals.first.insert(event_1_iterator+1,event_2_iterator,event_1_iterator+1);
		swapped_marginals.first.erase(event_2_iterator,event_1_iterator+1);
	}

	//Reassign the swapped array
	swapped_marginals.second = new_array_ptr;

}

/**
 * Utility to align marginals with events in the same order
 * Note that the implementation is probably not efficient but at least remains simple
 * Can align marginals with higher number of dimensions than the reference
 */
void align_marginal_array(const pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<double>>& reference_marginals , pair<list<pair<Rec_Event_name,size_t>>,shared_ptr<double>>& aligned_marginals){
	if(reference_marginals.first.size()>aligned_marginals.first.size()){
		throw runtime_error("Aligned marginals have less dimension than the reference marginals");
	}
	list<pair<Rec_Event_name,size_t>>::const_iterator reference_iterator = reference_marginals.first.begin();
	size_t counter = 0;
	while(reference_iterator != reference_marginals.first.end()){
		//Swap the reference event directly at the right position
		swap_events_order(reference_iterator->first , (aligned_marginals.first.begin()+counter)->first , aligned_marginals);
		++counter;
		++reference_iterator;
	}
}



