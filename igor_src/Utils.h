/*
 * Utils.h
 *
 *  Created on: Apr 9, 2015
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

#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <iostream>
#include "IntStr.h"
#include <memory>
#include <list>
#include <random>
#include <chrono>
#include<sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <unordered_map>


class Rec_Event;


enum Event_type {GeneChoice_t , Deletion_t , Insertion_t , Dinuclmarkov_t,Undefined_t};
enum Event_safety{VD_safe = 0  , DJ_safe = 1  , VJ_safe = 2 };
enum Seq_side{ Five_prime =0 , Three_prime = 1 , Undefined_side = 2 };
enum Seq_type {V_gene_seq = 0 , VD_ins_seq = 1 , D_gene_seq = 2 , DJ_ins_seq = 3 , J_gene_seq = 4 , VJ_ins_seq = 5};
enum Gene_class{V_gene=0 , VD_genes=1 , D_gene=2 , DJ_genes=3 , J_gene=4 , VJ_genes=5 , VDJ_genes=6 ,Undefined_gene=7 };
enum Fileformat{CSV_f,FASTA_f,TXT_f, FASTQ_f};
enum Int_nt { int_A = 0 , int_C = 1 , int_G = 2 , int_T = 3 , int_R = 4 , int_Y = 5 , int_K = 6 , int_M = 7 , int_S = 8 ,
				int_W = 9 , int_B = 10 , int_D = 11 , int_H = 12 , int_V = 13 , int_N = 14};


Gene_class str2GeneClass(const std::string);
std::string to_string(const Gene_class);
Seq_side str2SeqSide(const std::string);
std::string to_string(const Seq_side);

std::ostream& operator<<(std::ostream& , Gene_class);
std::ostream& operator<<(std::ostream& , Seq_side);
std::string operator+(const std::string& , Gene_class );
std::string operator+(const std::string& , Seq_side );
std::string operator+(const std::string& , Event_type );

//Type used to describe the array of doubles containing the marginals values
typedef std::unique_ptr<long double []> Marginal_array_p;

//Type used as key for unordered map since Rec_event cannot be instantiated
typedef std::string Rec_Event_name;

//Type used for offset of alignmed sequences in sequence_offsets maps. Used to characterize the beginning and the end of a sequence on the data sequence
typedef int Seq_Offset;
typedef Int_Str* Int_Str_ptr;

//Typedef used for getting the next event ptr
//typedef std::shared_ptr<Rec_Event> Next_event_ptr; //Does not work for some reason
typedef Rec_Event* Next_event_ptr;

/**
 * \brief Declare a null_delete function
 * \author Q.Marcou
 * This function is not performing any task, it's purpose is to supply a "null_delete" function
 * to prevent shared pointer objects created when passing Rec_Event or Error_rate objects pointers to model_parms
 * to be destroyed when the model_parms object is destroyed itself and the rec_event and error_rate objects
 * might still be of used (and if not prevent from a segfault error by trying to delete them twice)
 */
template<class T> struct null_delete{
	null_delete<T>(){}
	void operator()(T*&){};
};

/*
 * Declare a simple matrix class
 *
 */
template<typename T> struct Matrix {
public:
	Matrix():rows(0) , cols(0) , array_p(new T [0]){}
	Matrix(int m ,int n ): rows(m) , cols(n) , array_p(nullptr){
		if(m*n>0 and m>0){
			array_p = new T [m*n];
		}
	}
	Matrix(int m , int n , T arr[]):rows(m) , cols(n) , array_p(new T [m*n]){
		for(size_t i = 0 ; i != m*n ; i++){
			array_p[i] = arr[i];
		}
	}
	Matrix(int m , int n , std::vector<T> vect):rows(m) , cols(n) , array_p(new T [m*n]){
		for(size_t i = 0 ; i != m*n ; i++){
			array_p[i] = vect.at(i);
		}
	}
	Matrix(const Matrix<T>& other){
		//Provides deep copy of a matrix
		this->rows = other.rows;
		this->cols = other.cols;
		this->array_p = new T [rows*cols];
		for(int i = 0 ; i != rows*cols ; i++){
			this->array_p[i] = other.array_p[i];
		}
	}
	~Matrix(){
			delete [] array_p;
	}

	Matrix<T>& operator=(const Matrix& other){
		delete [] array_p;
		this->rows = other.rows;
		this->cols = other.cols;
		this->array_p = new T [rows*cols];
		for(int i = 0 ; i != rows*cols ; i++){
			this->array_p[i] = other.array_p[i];
		}
		return *this;
	}

	T& operator()(const int& i ,const int& j ){
		if( (i>rows-1) || (j>cols-1) ){
			throw std::length_error("Cannot access indices ["+std::to_string(i)+","+std::to_string(j)+"] with matrix dimensions ["+std::to_string(rows)+","+std::to_string(cols)+"]");
			std::cout<<"out_of range matrix coordinates: "<<rows<<"<"<<i<<" or "<<cols<<"<"<<j<<std::endl;
		}
		return array_p[i + rows*j];
	}

	const T& operator()(const int& i ,const int& j ) const{
		if( (i>rows-1) || (j>cols-1) ){
			throw std::length_error("Cannot access indices ["+std::to_string(i)+","+std::to_string(j)+"] with matrix dimensions ["+std::to_string(rows)+","+std::to_string(cols)+"]");
			std::cout<<"out_of range matrix coordinates: "<<rows<<"<"<<i<<" or "<<cols<<"<"<<j<<std::endl;
		}
		return array_p[i + rows*j];
	}

	T get_field(const int& i ,const int& j ) const{
		if( (i>rows-1) || (j>cols-1) ){
			throw std::length_error("Cannot access indices ["+std::to_string(i)+","+std::to_string(j)+"] with matrix dimensions ["+std::to_string(rows)+","+std::to_string(cols)+"]");
			std::cout<<"out_of range matrix coordinates: "<<rows<<"<"<<i<<" or "<<cols<<"<"<<j<<std::endl;
		}
		return array_p[i + rows*j];
	}

	//Accessors
	const int& get_n_rows() const {return rows;}
	const int& get_n_cols() const {return cols;}
private:
	int rows;
	int cols;
	T* array_p ;
};

template <typename T>
std::ostream& operator<<(std::ostream& stream ,const Matrix<T>& mat){
	stream<<mat.get_n_rows()<<"x"<<mat.get_n_cols()<<" Matrix"<<std::endl;
	for(int j = 0 ; j!= mat.get_n_cols() ; ++j){
		for(int i = 0 ; i!= mat.get_n_rows() ; ++i){
			if(i!=0){
				stream<<" ";
			}
			stream<<mat.get_field(i,j);
		}
		stream<<std::endl;
	}
	return stream;
}


/*
 * This class provides a fast alternative to unordered_map<Seq_type,string*> for the constructed_sequences objects
 * Change this and give some kind of matrix with memory levels
 * Create a 0 size at first?
 * Get rid of it in the deletions
 */
template<typename K, typename V> class  Enum_fast_memory_map {
public:
		Enum_fast_memory_map(int defined_range):max_layer(0) , range(defined_range){
			value_ptr_arr =  new V [range];
			memory_layer_ptr = new int [range];
			for(size_t i = 0 ; i!=range ; ++i){
				this->memory_layer_ptr[i] = -1;
			}
			/*for(size_t i = 0 ; i!=range ; i++){
				str_ptr_arr[i] = nullptr;
			}*/
		}
		virtual ~Enum_fast_memory_map(){
			delete [] value_ptr_arr;
			delete [] memory_layer_ptr;
		}

		//Accessors
		V& operator[](const K& key){
			if(key>range-1){throw std::out_of_range("Unknown key in Enum_fast_memory_map::operator(Seq_type)");}
			if(memory_layer_ptr[key]>-1){
				return (*(value_ptr_arr+ key + memory_layer_ptr[key]*range));
			}
			else{
				memory_layer_ptr[key] = 0;
				return (*(value_ptr_arr+key));
			}

		}


		V& at(const K& key){
			if(key>range-1){
				throw std::out_of_range("Unknown key in Enum_fast_memory_map::operator(K key)");
			}
			else{
				if(memory_layer_ptr[key]>-1){
					return (*(value_ptr_arr + key + memory_layer_ptr[key]*range));
				}
				else{
					throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_map::at(const K& key)");
				}
			}
		}


		V& at(const K& key ,int memory_layer){
			if(key>range-1){
				throw std::out_of_range("Unknown seq type in Enum_fast_memory_map::operator(Seq_type)");
			}
			else{
				if(memory_layer<=(memory_layer_ptr[key]+1)){
					memory_layer_ptr[key] = memory_layer;
					return (*(value_ptr_arr + key + memory_layer*range));
				}
				else{
					throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_map::at( const K& key, int memory_layer)");
				}
			}
		}

		const V& at(const K& key ,int memory_layer) const{
			if(key>range-1){
				throw std::out_of_range("Unknown seq type in Enum_fast_memory_map::operator(Seq_type)");
			}
			else{
				if(memory_layer<=(memory_layer_ptr[key]+1)){
					memory_layer_ptr[key] = memory_layer;
					return (*(value_ptr_arr + key + memory_layer*range));
				}
				else{
					throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_map::at( const K& key, int memory_layer)");
				}
			}
		}


		int get_current_memory_layer(const K& key){
			return memory_layer_ptr[key];
		}

		void get_all_current_memory_layer(int* memory_layers_recipient){
			for(size_t i =0 ; i!=range ; ++i){
				memory_layers_recipient[i] = memory_layer_ptr[i];
			}
		}

		bool exist(const K& key){
			return memory_layer_ptr[key]>-1;
		}


		void request_memory_layer(const K& key){
			/*std::cout<<key<<std::endl;
			std::cout<<memory_layer_ptr[key]<<std::endl;*/
			if(key>range-1){
				throw std::out_of_range("Unknown key in Enum_fast_memory_map::request_memory_layer()");
			}
			//Get current memory layer at this position
				if(memory_layer_ptr[key]<max_layer){
					++memory_layer_ptr[key];
				}
				else{
					++max_layer;
					V* new_value_ptr = new V [range*(max_layer+1)];
					for(size_t i = 0 ; i!=range ; ++i){
						for(size_t j=0 ; j!=(max_layer) ; ++j){
							(*(new_value_ptr + i+j*range)) = (*(value_ptr_arr + i+j*range));
						}
					}
					delete [] value_ptr_arr;
					value_ptr_arr = new_value_ptr;
					++memory_layer_ptr[key];
				}
			}


		//Setters
		void set_value(const K& key ,const V& value , int memory_layer){
			if(key>range-1){throw std::out_of_range("Unknown seq type in Seq_type_str_p_map::operator(Seq_type)");}
				//Cannot fill memory layer without filling the ones downstream
				if(memory_layer<=(memory_layer_ptr[key]+1)){
					(*(value_ptr_arr + key + memory_layer*range)) = value;
					//Setting a value at a given layer invalidate upper layers
					memory_layer_ptr[key] = memory_layer;
				}
				else{
					throw std::out_of_range("Trying to access incorrect memory layer in Enum_fast_memory_map::set_value()");
				}
		}

		void multiply_all(double& prod_operand , int* memory_adresses){
			for(size_t i = 0 ; i!=range ; ++i){
/*				std::cout<<i<<std::endl;
				std::cout<<(*(value_ptr_arr + i))<<std::endl;
				std::cout<<(*(memory_adresses + i))*range<<std::endl;*/
				prod_operand *= value_ptr_arr[ i +memory_adresses[i]*range];
			}
		}

		void reset(){
			for(size_t i = 0 ; i!=range ; ++i){
				if(memory_layer_ptr[i]>-1){
					memory_layer_ptr[i]=0;
				}
			}
		}

		void init_first_layer(V value){
			for(size_t i = 0 ; i!=range ; ++i){
				if(memory_layer_ptr[i]>-1){
					throw std::runtime_error("First memory layer already initialized for key " + std::to_string(i) + " in Enum_fast_memory_map::init_first_layer");
				}
				else{
					value_ptr_arr[i] = value;
					memory_layer_ptr[i] = 0;
				}
			}
		}



protected:
	V* value_ptr_arr;
	int* memory_layer_ptr;
	int max_layer;
	size_t range; //= 6; //number of outcomes in Seq_type

};



typedef Enum_fast_memory_map<Seq_type,Int_Str_ptr> Seq_type_str_p_map;

typedef Enum_fast_memory_map<Event_safety,bool> Safety_bool_map;

typedef Enum_fast_memory_map<Seq_type,std::vector<int>*> Mismatch_vectors_map;

typedef Enum_fast_memory_map<int,size_t> Index_map;

typedef Enum_fast_memory_map<Seq_type,double> Downstream_scenario_proba_bound_map;

/*
	template<> class Enum_fast_memory_map<Seq_type ,Str_ptr>{
	Enum_fast_memory_map():Enum_fast_memory_map<Seq_type,Str_ptr>(6){};
	size_t range = 6;
};
*/


/*template<>
class Seq_type_str_p_map : public Enum_fast_memory_map<Seq_type,Str_ptr>{

};*/


/*
 * This class provides a fast alternative to unordered_map<Seq_type,string*> for the constructed_sequences objects
 * Change this and give some kind of matrix with memory levels
 * Create a 0 size at first?
 * Get rid of it in the deletions
 */
template<typename K1, typename K2 , typename V> class  Enum_fast_memory_dual_key_map {
public:
		Enum_fast_memory_dual_key_map(size_t Key1_range , size_t Key2_range):max_layer(0) , range_key1(Key1_range) , range_key2(Key2_range){
			total_range = range_key1*range_key2;
			value_ptr_arr =  new V [total_range];
			memory_layer_ptr = new int [total_range];
			for(size_t i = 0 ; i!=total_range ; ++i){
				this->memory_layer_ptr[i] = -1;
			}
			/*for(size_t i = 0 ; i!=range ; i++){
				str_ptr_arr[i] = nullptr;
			}*/
		}
		virtual ~Enum_fast_memory_dual_key_map(){
			delete [] value_ptr_arr;
			delete [] memory_layer_ptr;
		}

		//Accessors
		//Cannot use [] with more than one argument
		/*V& operator[](const K1& key1 , const K2& key2){
			if(key1>range_key1-1){throw std::out_of_range("Unknown key1 in Enum_fast_memory_map::operator(Seq_type)");}
			if(key2>range_key2-1){throw std::out_of_range("Unknown key2 in Enum_fast_memory_map::operator(Seq_type)");}
			if(memory_layer_ptr[key1+range_key1*key2]>-1){
				return (*(value_ptr_arr+ key1+range_key1*key2 + memory_layer_ptr[key1+range_key1*key2]*total_range));
			}
			else{
				memory_layer_ptr[key1+range_key1*key2] = 0;
				return (*(value_ptr_arr+key1+range_key1*key2));
			}

		}*/


		V& at(const K1& key1 , const K2& key2){
			if(key1>range_key1-1){
				throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key_map::at()");
			}
			else if(key2>range_key2-1){
				throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
			}
			else{
				if(memory_layer_ptr[key1+range_key1*key2]>-1){
					return (*(value_ptr_arr + key1+range_key1*key2 + memory_layer_ptr[key1+range_key1*key2]*total_range));
				}
				else{
					throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
				}
			}
		}

		const V& at(const K1& key1 , const K2& key2) const{
			if(key1>range_key1-1){
				throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key_map::at()");
			}
			else if(key2>range_key2-1){
				throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
			}
			else{
				if(memory_layer_ptr[key1+range_key1*key2]>-1){
					return (*(value_ptr_arr + key1+range_key1*key2 + memory_layer_ptr[key1+range_key1*key2]*total_range));
				}
				else{
					throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
				}
			}
		}

		V& at(const K1& key1 , const K2& key2 , int memory_layer){
			if(key1>range_key1-1){
				throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key__map::at()");
			}
			else if(key2>range_key2-1){
				throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
			}
			else{
				if(memory_layer<=(memory_layer_ptr[key1+range_key1*key2]+1)){
					memory_layer_ptr[key1+range_key1*key2] = memory_layer;
					return (*(value_ptr_arr + key1+range_key1*key2 + memory_layer*total_range));
				}
				else{
					throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
				}
			}
		}

		const V& at(const K1& key1 , const K2& key2 , int memory_layer) const{
			if(key1>range_key1-1){
				throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key__map::at()");
			}
			else if(key2>range_key2-1){
				throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
			}
			else{
				if(memory_layer<=(memory_layer_ptr[key1+range_key1*key2]+1)){
					memory_layer_ptr[key1+range_key1*key2] = memory_layer;
					return (*(value_ptr_arr + key1+range_key1*key2 + memory_layer*total_range));
				}
				else{
					throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
				}
			}
		}


		int get_current_memory_layer(const K1& key1 , const K2& key2){
			return memory_layer_ptr[key1+range_key1*key2];
		}


		void request_memory_layer(const K1& key1 , const K2& key2){
			/*std::cout<<key<<std::endl;
			std::cout<<memory_layer_ptr[key]<<std::endl;*/
			if(key1>range_key1-1){
				throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key__map::request_memory_layer()");
			}
			if(key2>range_key2-1){
				throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::request_memory_layer()");
			}
			//Get current memory layer at this position
				if(memory_layer_ptr[key1+range_key1*key2]<max_layer){
					++memory_layer_ptr[key1+range_key1*key2];
				}
				else{
					++max_layer;
					V* new_value_ptr = new V [total_range*(max_layer+1)];
					for(size_t i = 0 ; i!=total_range ; ++i){
						for(size_t j=0 ; j!=(max_layer) ; ++j){
							(*(new_value_ptr + i+j*total_range)) = (*(value_ptr_arr + i+j*total_range));
						}
					}
					delete [] value_ptr_arr;
					value_ptr_arr = new_value_ptr;
					++memory_layer_ptr[key1+range_key1*key2];
				}
			}


		//Setters
		void set_value(const K1& key1 , const K2& key2 ,V value , int memory_layer){
			if(key1>range_key1-1){throw std::out_of_range("Unknown key1 in Seq_type_str_p_map::set_value()");}
			if(key2>range_key2-1){throw std::out_of_range("Unknown key2 in Seq_type_str_p_map::set_value()");}
				//Cannot fill memory layer without filling the ones downstream
				if(memory_layer<=(memory_layer_ptr[key1+range_key1*key2]+1)){
					(*(value_ptr_arr + key1+range_key1*key2 + memory_layer*total_range)) = value;
					//Setting a value at a given layer invalidate upper layers
					memory_layer_ptr[key1+range_key1*key2] = memory_layer;
				}
				else{
					throw std::out_of_range("Trying to access incorrect memory layer in Enum_fast_memory_dual_key__map::set_value()");
				}
		}


protected:
	V* value_ptr_arr;
	int* memory_layer_ptr;
	int max_layer;
	size_t range_key1; //= 6; //number of outcomes in Seq_type
	size_t range_key2;
	size_t total_range;

};

typedef Enum_fast_memory_dual_key_map<Seq_type,Seq_side,Seq_Offset> Seq_offsets_map;







/*
 * Defining a hash functions for Rec_Event, Gene_class and pair<Gene_class,Seq_side>
 */
 namespace std{
	 /*
 	 template<>
	 struct hash<Rec_Event>{
		inline std::size_t operator()(const Rec_Event& event) const{ //TODO inline?
			return  (((hash<int>()(event.get_class())
					^(hash<int>()(event.get_side())<<1 )) >>1)
					^(hash<int>()(event.get_priority())<<1)>>1)
					^(hash<int>()(event.get_realizations_map().size())<<1);
			//Note : only consider the size of the realization map and not what it contains for speed purposes
			//this should be enough to ensure no collisions
		}
	 };
	 */


/*
	 template<>
	 struct hash<Rec_Event*>{
		 std::size_t operator()(const Rec_Event*& event_point) const{
			 return hash<Rec_Event>()(*event_point);
		 }
	 };
	 */

	 template<>
	 	 struct hash<Seq_type>{
	 		 std::size_t operator()(const Seq_type& seq_t) const{
	 		 			return  hash<int>()(seq_t);
	 		 		}
	 	 };

	template<>
		 struct hash<Gene_class>{
			 std::size_t operator()(const Gene_class& gene) const{
						return  hash<int>()(gene);
					}
		 };

	 template<>
	 struct hash<std::pair<Gene_class,Seq_side>>{
		 std::size_t operator()(const pair<Gene_class,Seq_side>& gene_pair) const{
		 		 		return  (hash<Gene_class>()(gene_pair.first)
		 		 				^(hash<int>()(gene_pair.second) <<1))>>1;
		 		 	}
	 };

	 template<>
	 struct hash<std::tuple<Event_type,Gene_class,Seq_side>>{
		std::size_t operator()(const std::tuple<Event_type,Gene_class,Seq_side>& event_triplet) const{
			Event_type ev_type;
			Gene_class g_class;
			Seq_side s_side;
			std::tie(ev_type,g_class,s_side) = event_triplet;
			return ((hash<int>()(ev_type)
					^(hash<int>()(g_class) <<1)>>1)
					^(hash<int>()(s_side) <<1));

		}
	 };

	 template<>
	 struct hash<std::pair<Seq_type,Seq_side>>{
		std::size_t operator()(const std::pair<Seq_type,Seq_side> seq_pair) const{
			return  (hash<int>()(seq_pair.first)
					 		 				^(hash<int>()(seq_pair.second) <<1))>>1;
		}
	 };

	 template<>
	 struct hash<Event_safety>{
		 std::size_t operator()(const Event_safety ev_saf) const{
			 return (hash<int>()(ev_saf));
		 }
	 };
 }

 struct D_position_comparator{
 	 bool operator()(std::tuple<std::string,int,int,double> position_1 , std::tuple<std::string,int,int,double> position_2 ){
 		 return std::get<3>(position_1) > std::get<3>(position_2);
 	 }
 };

 struct inverse_offset_comparator{
	 bool operator()(const std::pair<std::shared_ptr<const Rec_Event>,int>& inv_offset_1 , const std::pair<std::shared_ptr<const Rec_Event>,int>& inv_offset_2){
		 return inv_offset_1.second < inv_offset_2.second;
	 }
 };

std::vector<std::string> extract_string_fields(const std::string,const std::string);

void show_progress_bar(std::ostream&,double, std::string prefix_message = "", size_t progress_bar_size = 70);
void close_progress_bar(std::ostream&, std::string prefix_message = "", size_t progress_bar_size = 70);
uint64_t draw_random_64bits_seed();


typedef std::unordered_map<std::string,std::string> UMCodonTable;



std::string translate(std::string seq);

#endif /* UTILS_H_ */
