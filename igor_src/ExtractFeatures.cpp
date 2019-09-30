/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GeneFeatures.cpp
 * Author: alfaceor
 * 
 * Created on 20 de Setembro de 2019, 12:41
 */

#include "ExtractFeatures.h"

ExtractFeatures::ExtractFeatures() {
}

ExtractFeatures::ExtractFeatures(const ExtractFeatures& orig) {
}

ExtractFeatures::~ExtractFeatures() {
}
/**
 * \brief Function to map the templateID and the sequence for V and J.
 * \param v_genomic vector with the genomic templates for V genes
 * \param j_genomic vector with the genomic templates for J genes
 * \return
 */
void ExtractFeatures::load_VJgenomicTemplates(vector<pair<string,string>> v_genomic, vector<pair<string,string>> j_genomic) {
	for (auto it = v_genomic.begin(); it != v_genomic.end(); ++it){
		UMap_v_genomic[(*it).first] = (*it).second;
	}
	for (auto it = j_genomic.begin(); it != j_genomic.end(); ++it){
		UMap_j_genomic[(*it).first] = (*it).second;
	}
	
}

/**
 * \brief load data files into GeneFeatures functor class
 * @param flnV_CDR3_anchors CDR3 anchors filename for V genes
 * @param flnJ_CDR3_anchors CDR3 anchors filename for J genes
 */
void ExtractFeatures::load_VJanchors(string flnV_CDR3_anchors, string flnJ_CDR3_anchors){
	UMap_v_CDR3_anchors = read_gene_anchors_csv(flnV_CDR3_anchors);
	UMap_j_CDR3_anchors = read_gene_anchors_csv(flnJ_CDR3_anchors);
}

//void ExtractFeatures::set_indexed_seqlist(vector<pair<const int,const string> >* pointer){
//	p_indexed_seqlist = pointer;
//}

void ExtractFeatures::set_sorted_alignments(unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>>* pointer){
	p_sorted_alignments = pointer;
}

/**
 * \brief Get the CDR3 as an instance of CDR3SeqData given the V and J Alignment_data structs.
 * \param seq_index
 * \param V_alignment
 * \param J_alignment
 * \return CDR3SeqData
 */
CDR3SeqData ExtractFeatures::extractCDR3(int seq_index){
	string seq_str = (*p_sorted_alignments)[seq_index].first;
	CDR3SeqData cdr3 = CDR3SeqData();
	cdr3.seq_index = seq_index;
	
	///////////////
	// V anchor  //
	///////////////
	// Get all the V alignments corresponding to seq_index
	vector<Alignment_data> Vec_V_Alignment_data = ((*p_sorted_alignments)[seq_index].second )[V_gene];
	// If V alignment was found then get the corresponding V CDR3 anchor
	if ( Vec_V_Alignment_data.size() > 0 ){
		Alignment_data v_alig = Vec_V_Alignment_data.front();
		cdr3.v_anchor  = getVAnchor4Seq(seq_str, v_alig);
	}
	
	///////////////
	// J anchor  //
	///////////////
	// Get all the J alignments corresponding to seq_index
	vector<Alignment_data> Vec_J_Alignment_data = ((*p_sorted_alignments)[seq_index].second )[J_gene];
	// If V alignment was found then get the corresponding V CDR3 anchor
	if ( Vec_J_Alignment_data.size() > 0 ){
		Alignment_data j_alig = Vec_J_Alignment_data.front();
		cdr3.j_anchor  = getJAnchor4Seq(seq_str, j_alig);
	}
	
	return cdr3;

}

/**
 * Calculate the V anchor in reference to the input sequence
 * \param seq_str
 * \param v_alig
 * \return cdr3_v_read_anch the anchor in reference to the input sequence.
 */
int ExtractFeatures::getVAnchor4Seq(string seq_str, Alignment_data v_alig){
	string v_gene_name    = "";
	string v_gene_str     = "";

	int    v_ins_size     = 0;
	int    v_dels_size    = 0;

	int cdr3_v_gene_anch;

	v_gene_name = v_alig.gene_name;
	v_gene_str  = UMap_v_genomic[v_alig.gene_name];
	v_ins_size  = distance(v_alig.insertions.begin(), v_alig.insertions.end() );
	v_dels_size = distance(v_alig.deletions .begin(), v_alig.deletions .end() );

	// Get the anchor from map and correct them.
	cdr3_v_gene_anch = UMap_v_CDR3_anchors[v_alig.gene_name];
	int v_ins_correction = std::count_if(
			v_alig.insertions.begin(), v_alig.insertions.end(), 
			// Lambda function for condition
			[cdr3_v_gene_anch](int inss){ 
						return (inss <= cdr3_v_gene_anch); 
			}
	);
	cdr3_v_gene_anch = cdr3_v_gene_anch + v_ins_correction;

	int cdr3_v_read_anch = cdr3_v_gene_anch + v_alig.offset;
	int dels_correction = std::count_if(
			v_alig.deletions.begin(), v_alig.deletions.end(), 
			// Lambda function for condition
			[cdr3_v_read_anch](int dels){ 
						return (dels <= cdr3_v_read_anch); 
			}
		);
	cdr3_v_read_anch = cdr3_v_read_anch + dels_correction; // ins_size before the cdr3 a
	
	return cdr3_v_read_anch;
}


/**
 * Calculate the J anchor in reference to the input sequence
 * \param seq_str
 * \param j_alig
 * \return cdr3_j_read_anch the anchor in reference to the input sequence.
 */
int ExtractFeatures::getJAnchor4Seq(string seq_str, Alignment_data j_alig){
	string j_gene_name    = "";
	string j_gene_str     = "";

	int    j_ins_size     = 0;
	int    j_dels_size    = 0;

	int cdr3_j_gene_anch;

	
	j_gene_name = j_alig.gene_name;
	j_gene_str  = UMap_j_genomic[j_alig.gene_name];
	j_ins_size  = distance(j_alig.insertions.begin(), j_alig.insertions.end() );
	j_dels_size = distance(j_alig.deletions .begin(), j_alig.deletions .end() );

	// Get the anchor from map and correct them.
	cdr3_j_gene_anch = UMap_j_CDR3_anchors[j_alig.gene_name];
	int j_ins_correction = std::count_if(
			j_alig.insertions.begin(), j_alig.insertions.end(), 
			// Lambda function for condition
			[cdr3_j_gene_anch](int inss){ 
						return (inss <= cdr3_j_gene_anch); 
			}
	);
	cdr3_j_gene_anch = cdr3_j_gene_anch + j_ins_correction;

	int cdr3_j_read_anch = cdr3_j_gene_anch + j_alig.offset;
	int j_dels_correction = std::count_if(
			j_alig.deletions.begin(), j_alig.deletions.end(), 
			// Lambda function for condition
			[cdr3_j_read_anch](int dels){ 
						return (dels <= cdr3_j_read_anch); 
			}
		);

	// In order to get the phelanine or triptophan in the sequence +3.
	cdr3_j_read_anch = cdr3_j_read_anch + j_dels_correction + 3; // ins_size before the cdr3 anchor
	return cdr3_j_read_anch;
}


/**
 * \brief Function to generate a CDR3 string line to be printed on a file.
 * \param cdr3InputSeq is CDR3SeqData instance 
 * \return csvline line string to be printed on a file.
 */
string ExtractFeatures::generateCDR3_csv_line(CDR3SeqData cdr3InputSeq){
	// seq_index;v_anchor;j_anchor;CDR3nt;CDR3aa
	string strCSVdelimiter = ";";
	string seq_str = (*p_sorted_alignments)[cdr3InputSeq.seq_index].first;
	
	stringstream sstm;
	sstm.str("");
	sstm << cdr3InputSeq.seq_index << strCSVdelimiter;
	
	
	bool bNoV = ((cdr3InputSeq.v_anchor < 0 ) or (cdr3InputSeq.v_anchor > seq_str.size() ) );
	bool bNoJ = ((cdr3InputSeq.j_anchor < 0 ) or (cdr3InputSeq.j_anchor > seq_str.size() ) );
	if ( bNoV or bNoJ ){
		sstm << strCSVdelimiter;
		sstm << strCSVdelimiter;
		sstm << strCSVdelimiter;
	}else{
		sstm << cdr3InputSeq.v_anchor  << strCSVdelimiter;
		sstm << cdr3InputSeq.j_anchor  << strCSVdelimiter;
		string strCDR3 = seq_str.substr(cdr3InputSeq.v_anchor, cdr3InputSeq.j_anchor - cdr3InputSeq.v_anchor);
		sstm << strCDR3 << strCSVdelimiter;
		sstm << translate(strCDR3);
	}
	
	return ( ""+ sstm.str() );
	
}

