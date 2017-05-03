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


unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_inverse_offset_map(const Model_Parms& model_parms){
	queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
	return get_inverse_offset_map(model_parms,model_queue);
}


unordered_map<Rec_Event_name,list<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_inverse_offset_map(const Model_Parms& model_parms, queue<shared_ptr<Rec_Event>> model_queue){
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
unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_offsets_map(const Model_Parms& model_parms){
	queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
	return get_offsets_map(model_parms,model_queue);
}

unordered_map<Rec_Event_name,vector<pair<shared_ptr<const Rec_Event>,int>>> Model_marginals::get_offsets_map(const Model_Parms& model_parms,queue<shared_ptr<Rec_Event>> model_queue){

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



unordered_map<Rec_Event_name,int> Model_marginals::get_index_map(const Model_Parms& model_parms){
	queue<shared_ptr<Rec_Event>> model_queue = model_parms.get_model_queue();
	return get_index_map( model_parms , model_queue);
}

unordered_map<Rec_Event_name,int> Model_marginals::get_index_map(const Model_Parms& model_parms , queue<shared_ptr<Rec_Event>> model_queue){


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
	const Event_realization& event_real = event_ptr->get_realizations_map().at(realization_name);
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
