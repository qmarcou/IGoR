/*
 * IntStr.h
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

#ifndef INTSTR_H_
#define INTSTR_H_

#include <vector>
#include <ostream>

class Int_Str: public std::vector<int> {

public:

	static const std::size_t npos = -1;

	Int_Str& operator+=(const Int_Str&);
	Int_Str& operator+=(const int&);
	Int_Str& operator+=( int&& );
	//Int_Str& operator+=(int);
	Int_Str& append(const Int_Str&);
	Int_Str& append(const int&);

	Int_Str operator+(const Int_Str&) const;
	Int_Str operator+(const int&) const;
	Int_Str operator+(int) const;

	Int_Str substr (std::size_t pos = 0, std::size_t len = npos ) const; //TODO correct this with aproper value

	using std::vector<int>::erase;
	Int_Str& erase(std::size_t pos , std::size_t len );



/*	Int_Str();
	//Int_Str (const std::string& str);//cannot be unambiguously identified
	Int_Str (const Int_Str&);
	//string (const string& str, size_t pos, size_t len = npos);//substring (3)
	//string (const char* s);//from c-string (4)
	//string (const char* s, size_t n);//from buffer (5)
	//string (size_t n, char c);//fill (6)
	//template <class InputIterator> string  (InputIterator first, InputIterator last);//range (7)
	//string (initializer_list<char> il);//initializer list (8)
	//string (string&& str) noexcept;//move (9)
	virtual ~Int_Str();
	Int_Str& operator=(const Int_Str&);*/

	//Iterators


/*
	//Capacity
	size_t max_size() const noexcept;
	//void resize();
	size_t capacity() const noexcept;
	void clear() noexcept;
	bool empty() const noexcept;
	size_t size() const noexcept;

	//Element access
	const int& operator[]() const;
	int& operator[]();
	const int& at() const;
	int& at() const;
	int& front();
	int& back();
*/


	//erase()

	//String operations
	//size_t find();
	//copy()
	//substr()
	//compare()

/*private:
	std::vector<int> int_vector;*/

};

namespace std{

template<>
	 struct hash<Int_Str>{
		std::size_t operator()(Int_Str const& int_str) const {
		  std::size_t seed = int_str.size();
		  for(auto& i : int_str) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		  }
		  return seed;
		}
	 };


}

#endif /* INTSTR_H_ */
