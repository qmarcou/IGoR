/*
 * Utils.cpp

 *
 *  Created on: Apr 9, 2015
 *      Author: quentin
 */

#include "Utils.h"

using namespace std;

/////Utilitaries
ostream& operator<<(ostream& os , Gene_class gc){
	switch(gc){
	case V_gene:os<<"V_gene";break;
	case VD_genes:os<<"VD_genes";break;
	case D_gene:os<<"D_gene";break;
	case DJ_genes:os<<"DJ_gene";break;
	case VJ_genes:os<<"VJ_gene";break;
	case J_gene:os<<"J_gene";break;
	case VDJ_genes:os<<"VDJ_genes";break;
	case Undefined_gene: os<<"Undefined_gene";break;

	default:
		throw invalid_argument("Unknown Gene_class in operator<< ");
	}
	return os;
}

ostream& operator<<(ostream& os , Seq_side ss){
	switch(ss){
	case Five_prime:os<<"Five_prime";break;
	case Three_prime:os<<"Three_prime";break;
	case Undefined_side:os<<"Undefined_side";break;
	}
	return os;
}

string operator+(const string& str , Gene_class gc){
	string next_str;
	switch(gc){
	case V_gene:next_str="V_gene";break;
	case VD_genes:next_str="VD_genes";break;
	case D_gene:next_str="D_gene";break;
	case DJ_genes:next_str="DJ_gene";break;
	case VJ_genes:next_str="VJ_gene";break;
	case J_gene:next_str="J_gene";break;
	case VDJ_genes:next_str="VDJ_genes";break;
	case Undefined_gene: next_str="Undefined_gene";break;
	}
	return str+next_str;
}

string operator+(const string& str , Seq_side ss){
	string next_str;
	switch(ss){
	case Five_prime:next_str="Five_prime";break;
	case Three_prime:next_str="Three_prime";break;
	case Undefined_side:next_str="Undefined_side";break;
	}
	return str+next_str;
}

string operator+(const string& str , Event_type et){
	string next_str;
	switch(et){
	case GeneChoice_t:next_str="GeneChoice";break;
	case Deletion_t:next_str="Deletion";break;
	case Insertion_t:next_str="Insertion";break;
	case Dinuclmarkov_t:next_str="DinucMarkov";break;
	}
	return str+next_str;
}

Gene_class str2GeneClass(string str){
	Gene_class gene_class;
	if(str == "V_gene"){gene_class = V_gene;}
	else if(str == "VD_genes"){gene_class = VD_genes;}
	else if(str == "D_gene"){gene_class = D_gene;}
	else if(str == "DJ_gene"){gene_class = DJ_genes;}
	else if(str == "VJ_gene"){gene_class = VJ_genes;}
	else if(str == "J_gene"){gene_class = J_gene;}
	else if(str == "VDJ_genes"){gene_class = VDJ_genes;}
	else if(str == "Undefined_gene"){gene_class = Undefined_gene;}
	else{
		throw runtime_error("Unknown Gene_class in str2GeneClass");
	}
	return gene_class;
}

Seq_side str2SeqSide(string str){
	Seq_side seq_side;
	if(str == "Five_prime"){seq_side = Five_prime;}
		else if(str == "Three_prime"){seq_side = Three_prime;}
		else if(str == "Undefined_side"){seq_side = Undefined_side;}
	else{
		throw runtime_error("Unknown event_side in str2SeqSide");
	}
	return seq_side;
}

//Seq_type_str_p_map::Seq_type_str_p_map(): Enum_fast_memory_map<Seq_type,Str_ptr>(6){}

//Implementing Seq_type_str_p_map methods

//template<typename K, typename V>

/*Enum_fast_memory_map::Enum_fast_memory_map(int defined_range) : max_layer(0),range(defined_range){


		str_ptr_arr =  new Str_ptr [range];
		memory_layer_ptr = new int [range];

		for(size_t i = 0 ; i!=range ; i++){
			str_ptr_arr[i] = nullptr;
		}

		for(size_t i = 0 ; i!=range ; i++){
			this->memory_layer_ptr[i] = -1;
		}
	};

Enum_fast_memory_map::~Enum_fast_memory_map(){
		delete [] str_ptr_arr;
		delete [] memory_layer_ptr;
	}


	 * Access operator without bound checks, once called the position is considered initialized

	Str_ptr& Seq_type_str_p_map::operator[](const Seq_type& seq_type){
		std::cout<<"operator[]"<<std::endl;
		for(size_t j=0 ; j!=(max_layer+1) ; j++){
						for(size_t i = 0 ; i!=range ; i++){
							std::cout<<(*(str_ptr_arr+i+range*j))<<";";
						}
						std::cout<<std::endl;
					}
					std::cout<<std::endl;


		if(seq_type>range-1){throw std::out_of_range("Unknown seq type in Seq_type_str_p_map::operator(Seq_type)");}
		if(memory_layer_ptr[seq_type]>-1){
			return (*(str_ptr_arr+seq_type + memory_layer_ptr[seq_type]*range));
		}
		else{
			memory_layer_ptr[seq_type] = 0;
			return (*(str_ptr_arr+seq_type));
		}

	}

	void update(){
		for(size_t i = 0 ; i!=range ; i++){
			if(tmp_modif_bool_array[i]){
				str_ptr_arr[i] = str_ptr_arr_tmp[i];
				tmp_modif_bool_array[i]=false;
				init_bool_array[i] = true;
			}
		}
	}


	 * This method returns the current memory layer at the given position

	int Seq_type_str_p_map::get_current_memory_layer(Seq_type seq_type){
		return memory_layer_ptr[seq_type];
	}


	 * Requesting a memory layer at a position adds a layer to that position

	void Seq_type_str_p_map::request_memory_layer(Seq_type seq_type){
		//Get current memory layer at this position
		if(memory_layer_ptr[seq_type]<max_layer){
			memory_layer_ptr[seq_type]++;
		}
		else{
			max_layer++;
			Str_ptr* new_str_ptr = new Str_ptr [range*(max_layer+1)];
			for(size_t i = 0 ; i!=range ; i++){
				for(size_t j=0 ; j!=(max_layer+1) ; j++){
					(*(new_str_ptr + i+j*range)) = (*(str_ptr_arr + i+j*range));
				}
			}
			delete [] str_ptr_arr;
			str_ptr_arr = new_str_ptr;
			memory_layer_ptr[seq_type]++;
		}
	}

	void Seq_type_str_p_map::set_value(Seq_type seq_type, Str_ptr value , int memory_layer){
		if(seq_type>range-1){throw std::out_of_range("Unknown seq type in Seq_type_str_p_map::operator(Seq_type)");}
		//Cannot fill memory layer without filling the ones downstream
		if(memory_layer<=(memory_layer_ptr[seq_type]+1)){
			(*(str_ptr_arr + seq_type + memory_layer*range)) = value;
			//Setting a value at a given layer invalidate upper layers
			memory_layer_ptr[seq_type] = memory_layer;
		}
		else{
			throw out_of_range("Trying to access incorrect memory layer");
		}

	}

	 * Access operator with bounds checks, prevent access to uninitialized positions

	Str_ptr& Seq_type_str_p_map::at(const Seq_type& seq_type){

		if(seq_type>range-1){
			throw std::out_of_range("Unknown seq type in Seq_type_str_p_map::operator(Seq_type)");
		}
		else{
			if(memory_layer_ptr[seq_type]>-1){
				return (*(str_ptr_arr + seq_type + memory_layer_ptr[seq_type]*range));
			}
			else{
				throw std::out_of_range("Trying to access unintialized position in Seq_type_str_p_map::at()");
			}
		}

	}


	Str_ptr& Seq_type_str_p_map::at(const Seq_type& seq_type , int memory_layer){

		if(seq_type>range-1){
			throw std::out_of_range("Unknown seq type in Seq_type_str_p_map::operator(Seq_type)");
		}
		else{
			if(memory_layer<=(memory_layer_ptr[seq_type]+1)){
				memory_layer_ptr[seq_type] = memory_layer;
				return (*(str_ptr_arr + seq_type + memory_layer*range));
			}
			else{
				throw std::out_of_range("Trying to access unintialized position in Seq_type_str_p_map::at()");
			}
		}

	}*/
/*
	 * Assignment operator,provides a shallow copy of the non tmp array

	Seq_type_str_p_map& operator=(const Seq_type_str_p_map& other){
		this->str_ptr_arr = other.str_ptr_arr;
		this->init_bool_array = other.init_bool_array;
		for(size_t i = 0 ; i!=range ; i++){
			if(other.tmp_modif_bool_array[i]){
				str_ptr_arr[i] = str_ptr_arr_tmp[i];

				init_bool_array[i] = true;
			}
		}
		return *this;
	}*/

	/*
	 * Copy constructor
	 * Provides a deep copy of the object

	Seq_type_str_p_map(const Seq_type_str_p_map& other){
		//Provides a deep copy of the "map"
		this->str_ptr_arr = new std::string* [range];

		this->init_bool_array = new bool [range];


		for(size_t i = 0 ; i!=range ; i++){
			this->str_ptr_arr[i] = other.str_ptr_arr[i];
			this->init_bool_array[i] = other.init_bool_array[i];
		}
	}*/
vector<string> extract_string_fields(string total_string ,string separator){
	vector<string> output_vector;
	int sep_index = -1;
	int next_sep_index = total_string.find(separator);
	while(next_sep_index>0){
		output_vector.emplace_back(total_string.substr(sep_index+1,next_sep_index));
		sep_index = next_sep_index;
		next_sep_index = total_string.find(separator,sep_index+1);
	}
	output_vector.emplace_back(total_string.substr(sep_index+1,string::npos));
	return output_vector;
}


