/*
 * IntStr.cpp
 *
 *  Created on: Jul 21, 2016
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

#include "IntStr.h"

using namespace std;

Int_Str& Int_Str::operator+=(const Int_Str& other){
	this->insert(this->end() , other.begin() , other.end());
	return *this;
}

Int_Str& Int_Str::operator +=(const int& a){
	this->push_back(a);
	return *this;
}

Int_Str& Int_Str::operator +=(int&& a ){
	this->push_back(a);
	return *this;
}

Int_Str& Int_Str::append(const Int_Str& other){
	return (*this)+= other;
}

Int_Str& Int_Str::append(const int& a){
	return (*this)+=a;
}

Int_Str Int_Str::operator+(const Int_Str& other) const{
	Int_Str new_int_str;
	new_int_str.reserve(this->size()+other.size());
	new_int_str+=(*this);
	new_int_str+=other;
	return new_int_str;
}

Int_Str Int_Str::substr(size_t pos /*= 0*/, size_t len /*= 99999999999999*/) const{
	Int_Str::const_iterator begin = this->begin() + pos;
	Int_Str::const_iterator last;
	if(len == npos){
		last = this->end();
	}
	else if( pos+len>=this->size()){
		last = this->end();
	}
	else{
		last = begin+len;
	}
	Int_Str new_int_str;
	new_int_str.assign(begin,last);
	return new_int_str;
}

Int_Str& Int_Str::erase(size_t pos , size_t len){
	Int_Str::iterator begin = this->begin() + pos;
	Int_Str::iterator last;
	if(len == npos){
		last = this->end();
	}
	else if( pos+len>=this->size()){
		last = this->end();
	}
	else{
		last = begin+len;
	}
	this->erase(begin,last);
	return *this;
}

ostream& operator<<(ostream& out , const Int_Str& int_str ){
	for(Int_Str::const_iterator iter = int_str.begin() ; iter!=int_str.end() ; ++iter){
		out<<(*iter);
	}
	return out;
}



/*
Int_Str::Int_Str(): int_vector() {
}

Int_Str::Int_Str(const Int_Str& other){
	this->int_vector = other.int_vector;
}

Int_Str::~Int_Str() {
	// TODO Auto-generated destructor stub
}

Int_Str& Int_Str::operator=(const Int_Str& other){
	this->int_vector = other.int_vector;
	return this;
}

size_t Int_Str::max_size() const noexcept{
	return this->int_vector.max_size();
}

size_t Int_Str::capacity() const noexcept{
	return this->int_vector.capacity();
}

void Int_Str::clear(){

}*/


