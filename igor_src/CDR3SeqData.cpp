/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CDR3SeqData.cpp
 * Author: alfaceor
 * 
 * Created on 7 de Junho de 2019, 15:46
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


