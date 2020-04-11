/*
 * File:   CDR3SeqData.cpp
 *
 *      Author: Carlos Olivares
 * 
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017- Quentin Marcou, 2019 - Carlos Olivares
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

#include <string>

#include "CDR3SeqData.h"

CDR3SeqData::CDR3SeqData() {
	seq_index = -1;
	v_anchor  = -1;
	j_anchor  = -1;
	CDR3nt    = "";
	CDR3aa    = "";
}

//CDR3SeqData::CDR3SeqData(int seq_index, int v_anchor, int j_anchor ) {
//	seq_index = seq_index;
//	v_anchor  = v_anchor;
//	j_anchor  = j_anchor;
//	CDR3nt    = "";
//	CDR3aa    = "";
//}

CDR3SeqData::CDR3SeqData(const CDR3SeqData& orig) {
}

CDR3SeqData::~CDR3SeqData() {
}


std::string CDR3SeqData::strData(){
	std::string delimiter = ";";
	return  std::to_string(seq_index) + delimiter +
					std::to_string(v_anchor ) + delimiter +
					std::to_string(j_anchor ) + delimiter +
					CDR3nt + delimiter +
					CDR3aa;	
}


