/*
 * IntStr.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: quentin
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

ostream& Int_Str::operator<<(ostream& out){
	for(vector<int>::const_iterator iter = this->begin() ; iter!=this->end() ; ++iter){
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


