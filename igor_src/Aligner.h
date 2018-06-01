/*
 * Aligner.h
 *
 *  Created on: Feb 16, 2015
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
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <forward_list>
#include <list>
#include <unordered_map>
#include <set>
#include <utility>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "Utils.h"
#include <omp.h>
#include <stdexcept>
#include <random>
#include <chrono>

#include "IntStr.h"


/**
 * \class Alignment_data Aligner.h
 * \brief Stores information on the alignment of one genomic template against the target
 * \author Q.Marcou
 * \version 1.0
 *
 * Stores information on the alignment of one genomic template against the target.
 * It contains:
 * - the gene name
 * - the offset of the alignment (index on the target sequence on which the first letter of the FULL genomic template aligns (can be negative or lie outside the target)
 * - 5' and 3' offset positions of the best alignment first and last aligned nucleotide
 * - insertions : indices on the TARGET of inserted nucleotides
 * - deletions : indices on the GENOMIC TEMPLATE of deleted nucleotides
 * - alignment length
 * - list of mismatches (that lie event outside the best alignment to allow IGoR to know mismatch positions in advance while exploring different deletions numbers)
 * - the alignment score
 */
struct Alignment_data {
	std::string gene_name;
	int offset;
	size_t five_p_offset;
	size_t three_p_offset;
	std::forward_list<int> insertions; //gap in the genomic sequence
	std::forward_list<int> deletions; //gap in the data sequence
	size_t align_length;
	mutable std::vector<int> mismatches;
	double score;

	Alignment_data(std::string gene , int off): gene_name(gene) , offset(off) , insertions(*(new std::forward_list<int>)) , deletions(*(new std::forward_list<int>)) , score(0) {}
	Alignment_data(int off, size_t five_p_off , size_t three_p_off , size_t align_len , std::forward_list<int> ins , std::forward_list<int> del , std::vector<int> mis , double alignment_score): gene_name(std::string()) , offset(off) , five_p_offset(five_p_off) , three_p_offset(three_p_off) , insertions(ins) , deletions(del) , align_length(align_len) , mismatches(mis) , score(alignment_score) {}
	Alignment_data(std::string gene , int off , size_t align_len , std::forward_list<int> ins , std::forward_list<int> del , std::vector<int> mis , double alignment_score): gene_name(gene) , offset(off) , insertions(ins) , deletions(del) , align_length(align_len) , mismatches(mis) , score(alignment_score) {}
	Alignment_data(std::string gene , int off, size_t five_p_off , size_t three_p_off , size_t align_len , std::forward_list<int> ins , std::forward_list<int> del , std::vector<int> mis , double alignment_score): gene_name(gene) , offset(off) , five_p_offset(five_p_off) , three_p_offset(three_p_off) , insertions(ins) , deletions(del) , align_length(align_len) , mismatches(mis) , score(alignment_score) {}

/*	bool operator<(const Alignment_data& align){
		//Hardcode to get the alignments in descending order using sort()
		return this->score > align.score;
	}*/

};

/**
 * \class Aligner Aligner.h
 * \brief A modified Smith-Waterman alignment class
 * \author Q.Marcou
 * \version 1.0
 *
 * The Aligner class allows to perform SW alignments according to the parameters (substitution matrix,gap penalty) supplied upon construction of the object.
 * The SW alignments matrix has been altered for V and J in order to allow for deletions on the deleted side only.
 * Alignments can be made in parallel using openMP
 *
 */
class Aligner {
public:
	Aligner();
	Aligner(Matrix<double>,int,Gene_class);
	virtual ~Aligner();

	// Single sequence alignments methods
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , int , int, bool=false);
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , int , int, std::set<std::string> , bool=false);
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , bool , int , int, bool=false);
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , bool , int , int, std::set<std::string> , bool=false);
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , std::unordered_map<std::string,std::pair<int,int>>, bool=false);
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , std::unordered_map<std::string,std::pair<int,int>>, std::set<std::string> , bool=false);
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , bool , std::unordered_map<std::string,std::pair<int,int>>, bool=false);
	std::forward_list<Alignment_data> align_seq(std::string , double , bool , bool , std::unordered_map<std::string,std::pair<int,int>> , std::set<std::string> , bool=false);

	// Multiple sequences alignments methods
	std::unordered_map<int,std::forward_list<Alignment_data>> align_seqs(std::vector<std::pair<const int , const std::string>> , double , bool);
	std::unordered_map<int,std::forward_list<Alignment_data>> align_seqs(std::vector<std::pair<const int , const std::string>> , double , bool , bool);
	std::unordered_map<int,std::forward_list<Alignment_data>> align_seqs(std::vector<std::pair<const int , const std::string>> , double , bool , int , int, bool=false);
	std::unordered_map<int,std::forward_list<Alignment_data>> align_seqs(std::vector<std::pair<const int , const std::string>> , double , bool , bool , int , int, bool=false);
	std::unordered_map<int,std::forward_list<Alignment_data>> align_seqs(std::vector<std::pair<const int , const std::string>> , double , bool , std::unordered_map<std::string,std::pair<int,int>>, bool=false);
	std::unordered_map<int,std::forward_list<Alignment_data>> align_seqs(std::vector<std::pair<const int , const std::string>> , double , bool , bool , std::unordered_map<std::string,std::pair<int,int>>, bool=false);
	void align_seqs( std::string , std::vector<std::pair<const int , const std::string>> , double , bool );
	void align_seqs( std::string , std::vector<std::pair<const int , const std::string>> , double , bool , bool );
	void align_seqs( std::string , std::vector<std::pair<const int , const std::string>> , double , bool , int , int, bool=false);
	void align_seqs( std::string , std::vector<std::pair<const int , const std::string>> , double , bool , bool , int , int, bool=false);
	void align_seqs( std::string , std::vector<std::pair<const int , const std::string>> , double , bool , std::unordered_map<std::string,std::pair<int,int>>, bool=false);
	void align_seqs( std::string , std::vector<std::pair<const int , const std::string>> , double , bool , bool , std::unordered_map<std::string,std::pair<int,int>>, bool=false);

	//I/O related methods
	void write_alignments_seq_csv(std::string , std::unordered_map<int,std::forward_list<Alignment_data>>);
	std::unordered_map<int,std::forward_list<Alignment_data>> read_alignments_seq_csv(std::string , double , bool);

	void set_genomic_sequences(std::vector< std::pair<std::string,std::string> >);
	int incorporate_in_dels( std::string& , std::string& , const std::forward_list<int> , const std::forward_list<int> , int );


private:
	std::forward_list<std::pair<std::string,std::string>> nt_genomic_sequences;
	std::forward_list<std::pair<std::string,Int_Str>> int_genomic_sequences;
	Matrix<double> substitution_matrix;
	int gap_penalty;
	Gene_class gene;
	bool local_align;
	bool flip_seqs;
	void sw_align_common(const Int_Str& ,const Int_Str& ,const int,const int , Matrix<double>& , Matrix<int>& , Matrix<int>& , Matrix<int>& , std::vector<int>& ,  std::vector<int>& , std::vector<int>&);
	std::list<std::pair<int,Alignment_data>> sw_align(const Int_Str& ,const Int_Str& , double , bool , int , int);
	std::unordered_map<std::string,std::pair<int,int>> build_genomic_bounds_map(int,int) const;

};


std::unordered_map<int,std::pair<std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>> read_alignments_seq_csv(std::string , Gene_class , double , bool , std::vector<std::pair<const int,const std::string>>);
std::unordered_map<int,std::pair<std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>> read_alignments_seq_csv(std::string , Gene_class , double , bool , std::vector<std::pair<const int,const std::string>>, std::unordered_map<int,std::pair<std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>>);
std::unordered_map<int,std::pair<std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>> read_alignments_seq_csv_score_range(std::string , Gene_class , double , bool , std::vector<std::pair<const int,const std::string>>);
std::unordered_map<int,std::pair<std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>> read_alignments_seq_csv_score_range(std::string , Gene_class , double , bool , std::vector<std::pair<const int,const std::string>>, std::unordered_map<int,std::pair<std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>>);
std::vector<std::tuple<int,std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>> map2vect (std::unordered_map<int,std::pair<std::string,std::unordered_map<Gene_class,std::vector<Alignment_data>>>>);
std::forward_list<std::pair<const int,const std::string>> read_indexed_seq_csv(std::string);
std::vector<std::pair<const int , const std::string>> read_indexed_csv(std::string);
std::vector<std::pair<const int,const std::string>> read_fasta(std::string);
std::vector<std::pair<std::string,std::string>> read_genomic_fasta(std::string);
std::vector<std::pair<const int,const std::string>> read_txt(std::string);
std::unordered_map<std::string,size_t> read_gene_anchors_csv(std::string,std::string separator= ";");
void write_indexed_seq_csv(std::string , std::vector<std::pair<const int,const std::string>>);
Int_Str nt2int(std::string);
bool comp_nt_int(const int& , const int&);
std::list<Int_nt> get_ambiguous_nt_list(const Int_nt&);
inline void write_single_seq_alignment( std::ofstream& , int , std::forward_list<Alignment_data> );
//Compare alignments (sort by score)
bool align_compare(Alignment_data , Alignment_data );
std::vector<std::pair<const int , const std::string>> sample_indexed_seq( std::vector<std::pair<const int , const std::string>>,const size_t);
Matrix<double> read_substitution_matrix(const std::string& , std::string sep=",");
std::tuple<bool,int,int> extract_min_max_genomic_templates_offsets(const std::unordered_map<std::string,std::pair<int,int>>& genomic_offset_bounds);
std::forward_list<Alignment_data> extract_best_gene_alignments(const std::forward_list<Alignment_data>&);

/*
	namespace substitution_matrices{
		//from: ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4
		static Matrix<int> nuc44_sub_matrix(4,4,{5,-4,-4,-4 , -4 ,5,-4,-4 , -4,-4,5,-4 , -4,-4,-4,5});


	}
	*/



#endif /* ALIGNER_H_ */
