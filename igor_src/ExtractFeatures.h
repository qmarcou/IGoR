

/* 
 * File:   ExtractFeatures.h
 * Author: alfaceor
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
    
//    void set_indexed_seqlist(vector<pair<const int, const string>>* pointer);
    void set_sorted_alignments(unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>>* pointer);
    
    
    CDR3SeqData generateCDR3(int seq_index);
    int getVAnchor4Seq(string seq_str, Alignment_data v_alig);
    int getJAnchor4Seq(string seq_str, Alignment_data j_alig);
    
    string generateCDR3_csv_line(CDR3SeqData cdr3InputSeq);
private:

};

#endif /* EXTRACTFEATURES_H */

