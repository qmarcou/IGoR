/*
 * File:   ExtractFeatures.h
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

#ifndef EXTRACTFEATURES_H
#define EXTRACTFEATURES_H
#include <unordered_map>
#include <vector>
#include "CDR3SeqData.h"
#include "Aligner.h"
#include <sstream>
#include "Utils.h"

using namespace std;

/**
 * \class ExtractFeatures ExtractFeatures.h
 * \brief Class to extract sequences features (e.g. CDR3) of sequences using alignment and V, J anchors information.
 * \author C. Olivares
 */
class ExtractFeatures {
public:
    ExtractFeatures();
    ExtractFeatures(const ExtractFeatures& orig);
    virtual ~ExtractFeatures();
    
    unordered_map<string, string> UMap_v_genomic;
    unordered_map<string, string> UMap_j_genomic;
    unordered_map<string, size_t> UMap_v_CDR3_anchors; 
    unordered_map<string, size_t> UMap_j_CDR3_anchors;
//    vector<pair<const int, const string>> *p_indexed_seqlist;
    unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>>  *p_sorted_alignments;
    
    void load_VJgenomicTemplates(vector<pair<string,string>> v_genomic, vector<pair<string,string>> j_genomic);
    void load_VJanchors(string flnV_CDR3_anchors, string flnJ_CDR3_anchors);
    void load_VJanchors(unordered_map<string, size_t>  flnV_CDR3_anchors, unordered_map<string, size_t>  flnJ_CDR3_anchors);
    
    void print_VgenomicTemplates();
    void print_JgenomicTemplates();
    
//    void set_indexed_seqlist(vector<pair<const int, const string>>* pointer);
    void set_sorted_alignments(unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>>* pointer);
    
    
    CDR3SeqData extractCDR3(int seq_index);
    int getVAnchor4Seq(string seq_str, Alignment_data v_alig);
    int getJAnchor4Seq(string seq_str, Alignment_data j_alig);
    
    string generateCDR3_csv_line(CDR3SeqData cdr3InputSeq);
private:

};

#endif /* EXTRACTFEATURES_H */

