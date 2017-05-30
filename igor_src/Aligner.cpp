/*
 * Aligner.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: quentin
 */

#include "Aligner.h"

using namespace std;

Aligner::Aligner() {
	// TODO Auto-generated constructor stub

}
/**
 * Constructor for the Aligner class
 * @sub_mat : substitution matrix
 * @gap_pen : sets the gap penalty (the gap penalty is linear)
 * @gene : Gene class of the gene aligned. V gene allows for deletions on the 3' side of the genomic template, J gene on the 5' , D gene and undefined allow deletion on both sides
 */
Aligner::Aligner(Matrix<double> sub_mat , int gap_pen , Gene_class gene): substitution_matrix(sub_mat) , gap_penalty(gap_pen){
	switch(gene){
	case V_gene:
		//Perform best alignment using all the right part of the genomic sequence
		local_align = false;
		flip_seqs = false;
		break;
	case D_gene:
		//Perform normal Smith Waterman alignment
		local_align = true;
		flip_seqs = false;
		break;
	case J_gene:
		//Perform best alignment using all the left part of the genomic sequence
		local_align = true;
		flip_seqs = false;
		break;
	case Undefined_gene:
		//Perform normal Smith Waterman alignment
		local_align = true;
		flip_seqs = false;
		break;
	default:
		throw runtime_error("Erroneous gene class for alignments");
		break;
	}

};

Aligner::~Aligner() {
	// TODO Auto-generated destructor stub
}
/*
 * This method reads sequences in a fasta file and return a vector of indexed sequences
 */
vector<pair<const int, const string>> read_fasta(string filename){
	//TODO Check for \r,\n\s stuff
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("File not found: "+filename);
	}
	string seq_str;
	string temp_str;
	int seq_count = -1;
	vector<pair<const int, const std::string>> sequence_vect;

	while (getline(infile,temp_str)){
		if(temp_str[temp_str.size()-1] == '\r'){
			temp_str.erase(temp_str.size()-1);
		}
		if(temp_str[0] == '>'){

			if(seq_count>(-1)){
				//Read sequences in upper case
				transform(seq_str.begin(),seq_str.end(),seq_str.begin(),::toupper);
				sequence_vect.push_back(pair<const int , const string> (seq_count , seq_str));
			}

			seq_str = string();
			seq_count++;
		}
		else{
			seq_str+=temp_str;
		}

	}
	if(seq_count>(-1)){
		//Read sequences in upper case
		transform(seq_str.begin(),seq_str.end(),seq_str.begin(),::toupper);
		sequence_vect.push_back(pair<const int , const string> (seq_count , seq_str));
	}

	return sequence_vect;
}

/*
 * This method reads genomic templates from fasta files
 */
vector<pair<string,string>> read_genomic_fasta(string filename){
	ifstream infile(filename);
		if(!infile){
			throw runtime_error("File not found: "+filename);
		}
		string seq_str;
		string seq_name_str;
		string temp_str = "";
		int seq_count = -1;
		vector<pair<string,std::string>> sequence_vect;

		while (getline(infile,temp_str)){
			if(temp_str[temp_str.size()-1] == '\r'){
				temp_str.erase(temp_str.size()-1);
			}
			if(temp_str[0] == '>'){

				if(seq_count>(-1)){
					//Convert all strings to upper case
					transform(seq_str.begin(),seq_str.end(),seq_str.begin(),::toupper);
					//Get rid of gaps inserted by IMGT to align the different genes
					string::iterator letter = seq_str.begin();
					while(letter != seq_str.end()){
						if( (*letter) =='.'){
							seq_str.erase(letter);
						}
						else if((*letter)=='\r'){
							//TODO Remove line return in sequences files in a proper way
							seq_str.erase(letter);
							break;
						}
						else{letter++;}
					}
					seq_name_str.erase(seq_name_str.begin()); //Get rid of the '>' at the beginning of the line
					sequence_vect.push_back(pair<string ,string> (seq_name_str , seq_str));
				}
				seq_name_str = temp_str;
				seq_str = string();
				seq_count++;
			}
			else{
				seq_str+=temp_str;
			}

		}
		if(seq_count>(-1)){
			//Convert all strings to upper case
			transform(seq_str.begin(),seq_str.end(),seq_str.begin(),::toupper);
			//Get rid of gaps inserted by IMGT to align the different genes
			string::iterator letter = seq_str.begin();
			while(letter != seq_str.end()){

				if( (*letter) =='.'){
					seq_str.erase(letter);
				}
				else if((*letter)=='\r'){
					//TODO Remove line return in sequences files in a proper way
					seq_str.erase(letter);
					break;
				}
				else{letter++;}
			}
			seq_name_str.erase(seq_name_str.begin());
			sequence_vect.push_back(pair<string ,string> (seq_name_str , seq_str));
		}

		return sequence_vect;
	}

/*
 * This method reads a file with one sequence per line and returns a vector of indexed sequences
 */
vector<pair<const int , const string>> read_txt(string filename){
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("File not found: "+filename);
	}
	int seq_count = 0;
	vector<pair<const int, const std::string>> sequence_vect;
	string seq_str;

	while(getline(infile,seq_str)){
		if(seq_str[seq_str.size()-1] == '\r'){
			seq_str.erase(seq_str.size()-1);
		}
		if(!seq_str.empty()){
			transform(seq_str.begin(),seq_str.end(),seq_str.begin(),::toupper);
			sequence_vect.push_back(pair<const int , const string >(seq_count , seq_str));
			seq_count++;
		}
	}

	return sequence_vect;
}
/*
 * This methods reads a file containing at each line the index of the sequence a semicolon and the actual sequence
 * @return : vector of indexed sequences vector<pair<const int , const string>>
 */
vector<pair<const int , const string>> read_indexed_csv(string filename){
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("File not found: "+filename);
	}
	string line_str;
	vector<pair<const int, const std::string>> sequence_vect;
	getline(infile,line_str);
	while(getline(infile,line_str)){
		size_t semi_col_index = line_str.find(";");
		int index = stoi(line_str.substr(0,semi_col_index));
		string seq_str = line_str.substr(semi_col_index+1 , string::npos);
		transform(seq_str.begin() , seq_str.end() , seq_str.begin() , ::toupper);
		sequence_vect.push_back(pair<const int , const string >(index , seq_str));
	}
	return sequence_vect;
}
/**
 * \brief This method aligns every genomic template to one sequence
 */
forward_list<Alignment_data> Aligner::align_seq(string nt_seq , double score_threshold , bool best_only , int min_offset , int max_offset){
	Int_Str int_seq = nt2int(nt_seq);
	forward_list<Alignment_data> alignment_list;// = *(new forward_list<Alignment_data>());
	for(forward_list<pair<string,Int_Str>>::const_iterator iter = int_genomic_sequences.begin() ; iter != int_genomic_sequences.end() ; iter++){
		list<pair<int,Alignment_data>> alignments = this->sw_align(int_seq , (*iter).second , score_threshold , best_only , min_offset , max_offset);
		//TODO quick and dirty fix for D genes alignments
		//alignment.second.gene_name = (*iter).first;
		//alignment_list.push_front(alignment.second);
		for(list<pair<int,Alignment_data>>::iterator jiter = alignments.begin() ; jiter != alignments.end() ; jiter++){
			(*jiter).second.gene_name = (*iter).first;
			alignment_list.push_front((*jiter).second);
		}
	}

	return alignment_list;
}

/*
 * Align sequences and hold them in memory
 */
unordered_map<int,forward_list<Alignment_data>> Aligner::align_seqs(vector<pair<const int , const string>> sequence_list , double score_threshold , bool best_only){
	unordered_map<int,forward_list<Alignment_data>> alignment_map = align_seqs(sequence_list , score_threshold , best_only , INT16_MIN , INT16_MAX);
	return alignment_map;
}

/*
 * Align sequences and hold them in memory
 */
unordered_map<int,forward_list<Alignment_data>> Aligner::align_seqs(vector<pair<const int , const string>> sequence_list , double score_threshold , bool best_only , int min_offset , int max_offset){
	unordered_map<int,forward_list<Alignment_data>> alignment_map; //= *(new unordered_map<int,forward_list<Alignment_data>>);

	int processed_seq_number = 0;

/*
 * Declaring parellel loop using OpenMP 4.0 standards
	#pragma omp declare reduction (merge:unordered_map<int,forward_list<Alignment_data>>:omp_out.insert(omp_in.begin(),omp_in.end()))
	#pragma omp parallel for schedule(dynamic) reduction(merge:alignment_map) shared(processed_seq_number)
*/

	//Declare parallel loop using OpenMP 3.1 standards
	#pragma omp parallel for schedule(dynamic) shared(processed_seq_number , alignment_map) //num_threads(1)
	for(vector<pair<const int , const string>>::const_iterator seq_it = sequence_list.begin() ; seq_it < sequence_list.end() ; seq_it++){
		forward_list<Alignment_data> seq_alignments = align_seq((*seq_it).second , score_threshold , best_only , min_offset , max_offset);
		#pragma omp critical(emplace_seq_alignments)
		{
			alignment_map.emplace((*seq_it).first , seq_alignments);
			//cout<<"Seq "<<processed_seq_number<<" processed"<<endl;
			processed_seq_number++;
		}
	}
	return alignment_map;
}

/*
 * Align sequences and write them on disk on the fly (avoids memory issues)
 */
void Aligner::align_seqs( string filename , vector<pair<const int , const string>> sequence_list , double score_threshold , bool best_only , int min_offset , int max_offset){
	unordered_map<int,forward_list<Alignment_data>> alignment_map; //= *(new unordered_map<int,forward_list<Alignment_data>>);

	string folder_path = filename.substr(0,filename.rfind("/")+1); //Get the file path
	ofstream align_infos_file(folder_path + "aligns_info.out",fstream::out | fstream::app); //Opens the file in append mode

	chrono::system_clock::time_point begin_time = chrono::system_clock::now();
	std::time_t tt;
	tt = chrono::system_clock::to_time_t ( begin_time );

	align_infos_file<<endl<<"================================================================"<<endl;
	align_infos_file<<"Alignments in file: "<<filename<<endl;
	align_infos_file<<"Date: "<< ctime(&tt)<<endl;
	align_infos_file<<"Score threshold = "<<score_threshold<<endl;
	align_infos_file<<"Best only = "<<best_only<<endl;
	align_infos_file<<"Min Offset = "<<min_offset<<endl;
	align_infos_file<<"Max Offset = "<<max_offset<<endl;
	align_infos_file<<"Gap penalty = "<<this->gap_penalty<<endl;
	align_infos_file<<"Substitution matrix:"<<endl;
	align_infos_file<<this->substitution_matrix<<endl;
	align_infos_file<<sequence_list.size()<<" sequences processed in ";

	ofstream outfile(filename);
	outfile<<"seq_index"<<";"<<"gene_name"<<";"<<"score"<<";"<<"offset"<<";"<<"insertions"<<";"<<"deletions"<<";"<<"mismatches"<<";"<<"length"<<";5_p_align_offset;3_p_align_offset"<<endl;

	int processed_seq_number = 0;

/*
 * Declaring parellel loop using OpenMP 4.0 standards
	#pragma omp declare reduction (merge:unordered_map<int,forward_list<Alignment_data>>:omp_out.insert(omp_in.begin(),omp_in.end()))
	#pragma omp parallel for schedule(dynamic) reduction(merge:alignment_map) shared(processed_seq_number)
*/

	//Declare parallel loop using OpenMP 3.1 standards
	#pragma omp parallel for schedule(dynamic) shared(processed_seq_number , alignment_map) //num_threads(1)
	for(vector<pair<const int , const string>>::const_iterator seq_it = sequence_list.begin() ; seq_it < sequence_list.end() ; seq_it++){
		try{
			forward_list<Alignment_data> seq_alignments = align_seq((*seq_it).second , score_threshold , best_only , min_offset , max_offset);

			#pragma omp critical(emplace_seq_alignments)
			{
				write_single_seq_alignment(outfile , (*seq_it).first , seq_alignments );
				//cout<<"Seq "<<processed_seq_number<<" processed"<<endl;
				++processed_seq_number;
			}
		}
		catch(exception& except){
			cout<<"Exception caught calling align_seq() on sequence:"<<endl;
			cout<<(*seq_it).first<<";"<<(*seq_it).second<<endl;
			cout<<endl;
			cout<<"Throwing exception now..."<<endl<<endl;
			cout<<except.what()<<endl;
			throw;
		}


	}

	chrono::duration<double> elapsed_time = chrono::system_clock::now() - begin_time;
	align_infos_file<<elapsed_time.count()<<" seconds"<<endl;

}

/*
 * Align sequences and write them on disk on the fly (avoids memory issues)
 */
void Aligner::align_seqs(string filename , vector<pair<const int , const string>> sequence_list , double score_threshold , bool best_only ){
	align_seqs(filename,sequence_list,score_threshold,best_only,INT16_MIN,INT16_MAX);
}

/*
 * Writes the indexed sequences as semicolon separated files with 2 fields:
 * @seq_index
 * @sequence
 */
void write_indexed_seq_csv(string filename , vector<pair<const int,const string>> indexed_seq_list){
	ofstream outfile(filename);
	outfile<<"seq_index"<<";"<<"sequence"<<endl;
	for(vector<pair<const int,const string>>::const_iterator iter = indexed_seq_list.begin() ; iter!=indexed_seq_list.end() ; iter++){
		outfile << (*iter).first<<";"<<(*iter).second<<endl;
	}
}
/*
 * Writes the alignment in a semicolon separated files with 5 fields:
 * @seq_index
 * @gene_name
 * @offset
 * @insertions: list of int coma separated surrounded by curly braces
 * @deletions: list of int coma separated surrounded by curly braces
 */
void Aligner::write_alignments_seq_csv(string filename , unordered_map<int,forward_list<Alignment_data>> indexed_alignments){
	ofstream outfile(filename);
	outfile<<"seq_index"<<";"<<"gene_name"<<";"<<"score"<<";"<<"offset"<<";"<<"insertions"<<";"<<"deletions"<<";"<<"mismatches"<<";"<<"length"<<endl;

	for(unordered_map<int,forward_list<Alignment_data>>::const_iterator iter = indexed_alignments.begin() ; iter != indexed_alignments.end() ; ++iter){
		write_single_seq_alignment(outfile , (*iter).first , (*iter).second);
	}
}

/*
 * This method writes the alignments for one sequence in the given stream
 */
void write_single_seq_alignment(ofstream& outfile , int seq_index , forward_list<Alignment_data> seq_alignments){
	for(forward_list<Alignment_data>::const_iterator jiter = seq_alignments.begin() ; jiter != seq_alignments.end() ; ++jiter){
		outfile<<seq_index<<";"<<(*jiter).gene_name<<";"<<(*jiter).score<<";"<<(*jiter).offset<<";{";
		for(forward_list<int>::const_iterator kiter = (*jiter).insertions.begin() ; kiter!=(*jiter).insertions.end() ; ++kiter){
			if(kiter==(*jiter).insertions.begin()){outfile<<(*kiter);}
			else{outfile<<","<<(*kiter);}
		}
		outfile<<"};{";
		for(forward_list<int>::const_iterator kiter = (*jiter).deletions.begin() ; kiter!=(*jiter).deletions.end() ; ++kiter){
			if(kiter==(*jiter).deletions.begin()){outfile<<(*kiter);}
			else{outfile<<","<<(*kiter);}
		}
		outfile<<"};{";//<<endl;
		for(vector<int>::const_iterator kiter = (*jiter).mismatches.begin() ; kiter!=(*jiter).mismatches.end() ; ++kiter){
			if(kiter==(*jiter).mismatches.begin()){outfile<<(*kiter);}
			else{outfile<<","<<(*kiter);}
		}
		outfile<<"};"<<(*jiter).align_length<<";"<<(*jiter).five_p_offset<<";"<<(*jiter).three_p_offset<<endl;
	}
}



/*
 * This method reads the indexed sequences from a given file(@filename)
 * The structure of the file is assumed to be the same as the one created by the Aligner::write_indexed_seq_csv method
 */
forward_list<pair<const int , const string>> read_indexed_seq_csv(string filename){
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("File not found: "+filename);
	}
	string line_str;
	forward_list<pair<const int , const string>> indexed_seq_list;
	//get rid of the first line
	getline(infile,line_str);
	while(getline(infile,line_str)){
		size_t scolon_index = line_str.find(';');
		indexed_seq_list.push_front(pair<const int,const string>(stoi(line_str.substr(0,scolon_index)),line_str.erase( 0 , (scolon_index+1) )));
	}
	return indexed_seq_list;
}



/*
 * This method reads the alignment data from a given file (@filename).
 * The structure of the file is assumed to be the same as the one created by the Aligner::write_alignments_seq_csv method
 */
unordered_map<int,vector<Alignment_data>> read_alignments_seq_csv(string filename , double score_threshold , bool allow_in_dels ){
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("File not found: "+filename);
	}
	string line_str;
	unordered_map<int,vector<Alignment_data>> indexed_alignments;
	//get rid of the first line
	getline(infile,line_str);
	while(getline(infile,line_str)){

		//find the semicolons in the line
		size_t index_sep = line_str.find(';');
		size_t name_sep = line_str.find(';',index_sep+1);
		size_t score_sep = line_str.find(';',name_sep+1);
		size_t off_sep = line_str.find(';',score_sep+1);
		size_t ins_sep = line_str.find(';',off_sep+1);
		size_t del_sep = line_str.find(';',ins_sep+1);
		size_t mism_sep = line_str.find(';',del_sep+1);

		int index = stoi(line_str.substr(0,index_sep));
		string gene_name = line_str.substr( (index_sep+1) , (name_sep-index_sep-1) );
		double score = stod(line_str.substr((name_sep+1) , (score_sep - name_sep -1)));
		int offset = stoi(line_str.substr( (score_sep+1) , (off_sep-score_sep-1) ));
		forward_list<int> insertions;
		forward_list<int> deletions;
		vector<int> mismatches;//TODO preallocate memory given the length of the string

		if(score < score_threshold){
			continue;
		}

		//Define a scope
		{
			//Get the index of insertions from comma separated integers surrounded by curly braces

			string ins_substr = line_str.substr( (off_sep+2) , (ins_sep - off_sep -3));//get rid of curly braces at the same time

			size_t comma_index =  ins_substr.find(',');
			if(comma_index!=string::npos){
				insertions.push_front(stoi(ins_substr.substr(0,(comma_index))));
				while(comma_index!=string::npos){
					size_t next_comma_index = ins_substr.find(',', (comma_index+1) );
					insertions.push_front(stoi(ins_substr.substr( (comma_index+1) , (next_comma_index - comma_index -1) )));
					comma_index = next_comma_index;

				}
			}
			else{
				if(!ins_substr.empty()){
					insertions.push_front(stoi(ins_substr));
				}
			}
		}

		{
			//Same with deletions

					string del_substr = line_str.substr( (ins_sep+2) , (del_sep - ins_sep -3));//get rid of curly braces at the same time

					size_t comma_index =  del_substr.find(',');
					if(comma_index!=string::npos){
						deletions.push_front(stoi(del_substr.substr(0,(comma_index))));
						while(comma_index!=string::npos){
							size_t next_comma_index = del_substr.find(',', (comma_index+1) );
							deletions.push_front(stoi(del_substr.substr( (comma_index+1) , (next_comma_index - comma_index -1) )));
							comma_index = next_comma_index;

						}
					}
					else{
						if(!del_substr.empty()){
							try{
								deletions.push_front(stoi(del_substr));
							}
							catch(exception& except){
								cout<<del_substr<<" cannot be casted as an integer in line:"<<endl;
								cout<<line_str<<endl;
								cout<<"Throwing exception now"<<endl;
								throw except;
							}
						}
					}

			}
		if(!allow_in_dels & (!insertions.empty() | !deletions.empty())){
			continue;
		}

		{
					//Same with mismatches
					string mismatch_substr;
					if(mism_sep==string::npos){
						//TODO remove this, this ensure compatibility with previous versions
						mismatch_substr = line_str.substr( (del_sep+2) , (line_str.size() - del_sep -3));//get rid of curly braces at the same time
					}
					else{
						mismatch_substr = line_str.substr( (del_sep+2) , (mism_sep - del_sep -3));//get rid of curly braces at the same time
					}


					size_t comma_index =  mismatch_substr.find(',');
					if(comma_index!=string::npos){
						mismatches.push_back(stoi(mismatch_substr.substr(0,(comma_index))));
						while(comma_index!=string::npos){
							size_t next_comma_index = mismatch_substr.find(',', (comma_index+1) );
							mismatches.push_back(stoi(mismatch_substr.substr( (comma_index+1) , (next_comma_index - comma_index -1) )));
							comma_index = next_comma_index;

						}
					}
					else{
						if(!mismatch_substr.empty()){
							mismatches.push_back(stoi(mismatch_substr));
						}
					}

					//FIXME read alignment length

			}

		indexed_alignments[index].push_back( Alignment_data(gene_name , offset , INT16_MIN , insertions , deletions , mismatches , score));

	}
	return indexed_alignments;
}

unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> read_alignments_seq_csv(string filename , Gene_class aligned_gene , double score_threshold , bool allow_in_dels , vector<pair<const int , const string>> indexed_sequences ){
	unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments;
	sorted_alignments = read_alignments_seq_csv(filename , aligned_gene , score_threshold , allow_in_dels , indexed_sequences , sorted_alignments);
	return sorted_alignments;
}

unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> read_alignments_seq_csv(string filename , Gene_class aligned_gene , double score_threshold , bool allow_in_dels , vector<pair<const int , const string>> indexed_sequences , unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments){
	unordered_map<int,vector<Alignment_data>> alignments = read_alignments_seq_csv(filename , score_threshold , allow_in_dels);
	for(vector<pair<const int , const string>>::const_iterator seq_it = indexed_sequences.begin() ; seq_it != indexed_sequences.end() ; ++seq_it){
		sorted_alignments[(*seq_it).first].second[aligned_gene] = alignments[(*seq_it).first];
		sorted_alignments[(*seq_it).first].first = (*seq_it).second;
	}
	//sort(alignments.begin() , alignments.end() , align_compare);
	return sorted_alignments;
}

unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> read_alignments_seq_csv_score_range(string filename , Gene_class aligned_gene , double score_range , bool allow_in_dels , vector<pair<const int , const string>> indexed_sequences ){
	unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments;
	sorted_alignments = read_alignments_seq_csv_score_range(filename , aligned_gene , score_range , allow_in_dels , indexed_sequences , sorted_alignments);
	return sorted_alignments;
}

unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> read_alignments_seq_csv_score_range(string filename , Gene_class aligned_gene , double score_range , bool allow_in_dels , vector<pair<const int , const string>> indexed_sequences , unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments){
	unordered_map<int,vector<Alignment_data>> alignments = read_alignments_seq_csv(filename , 0 , allow_in_dels);
	for(vector<pair<const int , const string>>::const_iterator seq_it = indexed_sequences.begin() ; seq_it != indexed_sequences.end() ; ++seq_it){
			vector<Alignment_data>& seq_alignments = alignments[(*seq_it).first];
			double max_score = -1;
			for(vector<Alignment_data>::const_iterator align_it = seq_alignments.begin() ; align_it != seq_alignments.end() ; ++align_it){
				if((*align_it).score>max_score){max_score=(*align_it).score;}
			}
			for(vector<Alignment_data>::iterator align_it = (seq_alignments.end()-1) ; align_it != (seq_alignments.begin() -1);--align_it){
				if((*align_it).score<(max_score-score_range)){seq_alignments.erase(align_it);}
			}
			sort(seq_alignments.begin() , seq_alignments.end() , align_compare);
			sorted_alignments[(*seq_it).first].second[aligned_gene] = alignments[(*seq_it).first];
			sorted_alignments[(*seq_it).first].first = (*seq_it).second;
	}
	return sorted_alignments;
}

vector<tuple<int,string,unordered_map<Gene_class,vector<Alignment_data>>>> map2vect(unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> alignments_map){
	vector<tuple<int,string,unordered_map<Gene_class,vector<Alignment_data>>>> alignmets_vect;
	for(unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>>::const_iterator seq_it = alignments_map.begin() ; seq_it != alignments_map.end() ; ++seq_it){
		alignmets_vect.emplace_back((*seq_it).first,(*seq_it).second.first , (*seq_it).second.second);
	}
	return alignmets_vect;
}

void Aligner::set_genomic_sequences(vector< pair <string,string> > nt_genomic_seq){
	this->nt_genomic_sequences = *(new forward_list<pair<string,string>>);
	this->int_genomic_sequences = *(new forward_list<pair<string,Int_Str>>);
	for(vector<pair<string,string>>::const_iterator iter = nt_genomic_seq.begin() ; iter != nt_genomic_seq.end() ; ++iter){
		nt_genomic_sequences.emplace_front((*iter).first,(*iter).second);
		int_genomic_sequences.emplace_front((*iter).first , nt2int((*iter).second));
	}
}

/*
 * This method will incorporate gaps('-') at the places where insertion or deletion occured both in the data and genomic sequence
 * A deletion will correspond to a gap introduced in the data sequence
 * An insertion to a gap in the genomic sequence
 * <\code>prev_dels<code> will indicate the shift in the data_seq offset induced by the introduction of previous deletions
 *
 * The method returns the shift induced by introducing the deletions of this alignment
 */
int Aligner::incorporate_in_dels(string& data_seq , string& genomic_seq , const forward_list<int> , const forward_list<int> , int prev_dels){

	return prev_dels;
}

/*
 * Convert nucleotide alphabet sequence into int sequence
 * Conventions are the same as nt2int matlab function
 */
Int_Str nt2int(string nt_sequence){
	Int_Str int_seq;
	for(size_t i=0 ; i!= nt_sequence.size() ; ++i){
		if(nt_sequence[i]=='A'){int_seq.append(int_A);}
		else if(nt_sequence[i]=='C'){int_seq.append(int_C);}
		else if(nt_sequence[i]=='G'){int_seq.append(int_G);}
		else if((nt_sequence[i]=='T') or (nt_sequence[i]=='U')){int_seq.append(int_T);}
		else if(nt_sequence[i]=='R'){int_seq.append(int_R);}
		else if(nt_sequence[i]=='Y'){int_seq.append(int_Y);}
		else if(nt_sequence[i]=='K'){int_seq.append(int_K);}
		else if(nt_sequence[i]=='M'){int_seq.append(int_M);}
		else if(nt_sequence[i]=='S'){int_seq.append(int_S);}
		else if(nt_sequence[i]=='W'){int_seq.append(int_W);}
		else if(nt_sequence[i]=='B'){int_seq.append(int_B);}
		else if(nt_sequence[i]=='D'){int_seq.append(int_D);}
		else if(nt_sequence[i]=='H'){int_seq.append(int_H);}
		else if(nt_sequence[i]=='V'){int_seq.append(int_V);}
		else if(nt_sequence[i]=='N'){int_seq.append(int_N);}
		else{
			cout<<"print:"<<nt_sequence<<endl;
			cout<<i<<endl;
			throw runtime_error("Unknown nucleotide: "+to_string(nt_sequence[i]) + "in string " + nt_sequence + "in Aligner::nt2int");
		}
	}
	return int_seq;
}

/**
 * This function compares nucleotides and output a boolean if they do not necessarily imply an error (ambiguous nucleotides are thus treated in a loose sense).
 */
bool inline comp_nt_int(const int& nt_1 , const int& nt_2){
	if(nt_1 != nt_2){
		if( (nt_1<4) & (nt_2<4)){
			return false;
		}
		else{
			switch(nt_1){
					case int_A:
						switch(nt_2){
							case int_R: case int_W: case int_M: case int_D: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_C:
						switch(nt_2){
							case int_Y: case int_S: case int_M: case int_B: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_G:
						switch(nt_2){
							case int_R: case int_S: case int_K: case int_B: case int_D: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_T:
						switch(nt_2){
							case int_Y: case int_W: case int_K: case int_B: case int_D: case int_H: case int_N:
								return true;
								break;
						}
						break;
					case int_R:
						switch(nt_2){
							case int_A: case int_G:
							case int_S: case int_W: case int_K: case int_M:
							case int_B: case int_D: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_Y:
						switch(nt_2){
							case int_C: case int_T:
							case int_S: case int_W: case int_K: case int_M:
							case int_B: case int_D: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_K:
						switch(nt_2){
							case int_G: case int_T:
							case int_R: case int_Y: case int_S: case int_W:
							case int_B: case int_D: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_M:
						switch(nt_2){
							case int_A: case int_C:
							case int_R: case int_Y: case int_S: case int_W:
							case int_B: case int_D: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_S:
						switch(nt_2){
							case int_G: case int_C:
							case int_R: case int_Y: case int_K: case int_M:
							case int_B: case int_D: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_W:
						switch(nt_2){
							case int_A: case int_T:
							case int_R: case int_Y: case int_K: case int_M:
							case int_B: case int_D: case int_H: case int_V: case int_N:
								return true;
								break;
						}
						break;
					case int_B:
						if(nt_2 != int_A){ //Of course this is in the hope that nt_2 is in the correct range of int
							return true;
						}
						break;
					case int_D:
						if(nt_2 != int_C){ //Of course this is in the hope that nt_2 is in the correct range of int
							return true;
						}
						break;
					case int_H:
						if(nt_2 != int_G){ //Of course this is in the hope that nt_2 is in the correct range of int
							return true;
						}
						break;
					case int_V:
						if(nt_2 != int_T){ //Of course this is in the hope that nt_2 is in the correct range of int
							return true;
						}
						break;
					case int_N:
						return true;
						break;
					default:
						throw runtime_error("Unknown nucleotide index: "+to_string(nt_1) + "in comp_nt_int()");

			}
			return false;
		}
	}
	else{
		return true;
	}
}

list<Int_nt> get_ambiguous_nt_list(const Int_nt& ambiguous_nt){
	list<Int_nt> nt_list;

	if ( (ambiguous_nt == int_A) or (ambiguous_nt == int_C) or  (ambiguous_nt ==int_G) or (ambiguous_nt == int_T)){
		nt_list.emplace_back(ambiguous_nt);
	}
	else{
		bool any_true=false;
		//Add an A for all cases implying an A	(do not add a break to allow for execution of other cases and possibly add more letters to the list)
		if((ambiguous_nt == int_R) or (ambiguous_nt == int_W) or (ambiguous_nt == int_M) or (ambiguous_nt == int_D)
				or (ambiguous_nt == int_H) or (ambiguous_nt ==int_V) or (ambiguous_nt ==int_N)){
			nt_list.emplace_back(int_A);
			any_true=true;
		}
		//Same for C
		if((ambiguous_nt == int_Y) or (ambiguous_nt == int_S) or (ambiguous_nt ==int_M) or (ambiguous_nt == int_B)
				or (ambiguous_nt ==int_H) or (ambiguous_nt == int_V) or (ambiguous_nt ==int_N)){
			nt_list.emplace_back(int_C);
			any_true=true;
		}
		//Same for G
		if((ambiguous_nt == int_R) or (ambiguous_nt ==int_S) or (ambiguous_nt ==int_K) or (ambiguous_nt ==int_B)
				or (ambiguous_nt ==int_D) or (ambiguous_nt ==int_V) or (ambiguous_nt ==int_N)){
			nt_list.emplace_back(int_G);
			any_true=true;
		}
		//Same for T
		if((ambiguous_nt == int_Y) or (ambiguous_nt ==int_W) or (ambiguous_nt ==int_K) or (ambiguous_nt ==int_B)
				or (ambiguous_nt ==int_D) or (ambiguous_nt ==int_H) or (ambiguous_nt ==int_N)){
			nt_list.emplace_back(int_T);
			any_true=true;
		}

		if( not any_true){
			throw runtime_error("Unknown nucleotide index: "+to_string(ambiguous_nt) + "in get_ambiguous_nt_list()");
		}

	}
	return nt_list;
}

/*
 * Randomly sample N indexed sequences from a given vector of indexed sequences
 */
vector<pair<const int , const string>> sample_indexed_seq(vector<pair<const int , const string>> indexed_seqs , const size_t sample_size){

	//Return an error if trying to sample more than the number of available sequences
	if(sample_size>indexed_seqs.size()){
		throw std::runtime_error("Trying to sample " + to_string(sample_size) + " sequences in a pool of " + to_string(indexed_seqs.size()) + " sequences in sample_indexed_seq()");
	}

	//Create seed for random generator
	//create a seed from timer
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point time = myclock::now();
	myclock::duration dur = myclock::time_point::max() - time;

	unsigned time_seed = dur.count();
	//Instantiate random number generator
	default_random_engine generator =  default_random_engine(time_seed);

	//Need to make a copy because of the constness in indexed_seqs
	vector<pair< int , string>> indexed_seqs_copy (indexed_seqs.begin(),indexed_seqs.end());

	shuffle(indexed_seqs_copy.begin(),indexed_seqs_copy.end(),generator);
	return vector<pair<const int , const string>>(indexed_seqs_copy.begin(),indexed_seqs_copy.begin()+sample_size);

}

void Aligner::sw_align_common(const Int_Str& int_data_sequence, const Int_Str& int_genomic_sequence ,const int i ,const int j , Matrix<double>& score_matrix , Matrix<int>& row_memory_matrix , Matrix<int>& col_memory_matrix , Matrix<int>& alignment_numb_tracker ,vector<int>& max_score, vector<int>& max_row_coord , vector<int>& max_col_coord){
	int genomic_gap_score = score_matrix(i,j-1) - gap_penalty;
	int data_gap_score = score_matrix(i-1,j) - gap_penalty;
	int subs_score = score_matrix(i-1,j-1) + substitution_matrix(int_data_sequence.at(i) , int_genomic_sequence.at(j));

	if ( (subs_score>= data_gap_score) & (subs_score>=genomic_gap_score) & ( (subs_score>0) | (!local_align) ) ){
		score_matrix(i,j) = subs_score;
		row_memory_matrix(i,j) = 1;
		col_memory_matrix(i,j) = 1;
		if(alignment_numb_tracker(i-1,j-1) == -1){
			alignment_numb_tracker(i,j) = max_score.size();
			max_score.push_back(subs_score);
			max_row_coord.push_back(i);
			max_col_coord.push_back(j);
		}
		else{
			alignment_numb_tracker(i,j) = alignment_numb_tracker(i-1,j-1);
		}

	}
	else if ( (data_gap_score>=genomic_gap_score) & ( (data_gap_score>0) | (!local_align) ) ){
		//This is in favor of deletions in the sequenced read instead of insertion.
		//No real rational behind it just to avoid undefined behavior of the aligner
		score_matrix(i,j) = data_gap_score;
		row_memory_matrix(i,j) = 1;
		col_memory_matrix(i,j) = 0;
		alignment_numb_tracker(i,j) = alignment_numb_tracker(i-1,j);
	}
	else if ( (genomic_gap_score>0) | (!local_align) ){
		score_matrix(i,j) = genomic_gap_score;
		row_memory_matrix(i,j) = 0;
		col_memory_matrix(i,j) = 1;
		alignment_numb_tracker(i,j) = alignment_numb_tracker(i,j-1);

	}
	else{
		score_matrix(i,j) = 0;
		row_memory_matrix(i,j) = 0;
		col_memory_matrix(i,j) = 0;
		//TODO check this
		alignment_numb_tracker(i,j) = alignment_numb_tracker(i-1,j-1);

	}
	//Keep max score in memory

	if(alignment_numb_tracker(i,j) != -1){
		if(score_matrix(i,j)>max_score[alignment_numb_tracker(i,j)]){
			int& tmp = alignment_numb_tracker(i,j);
			max_score[tmp] = score_matrix(i,j);
			max_row_coord[tmp] = i;
			max_col_coord[tmp] = j;
		}
	}


}
/*
 * Performs Smith-Waterman alignment between two sequences (translated to int sequence as a prior)
 * Output:
 * int:Alignment score
 * Alignment_data: comprises offset, insertions and deletions locations.
 * Note: the gene_name field of the Alignment_data object is left blank and should be completed in a higher level method
 */
list<pair<int,Alignment_data>> Aligner::sw_align(const Int_Str& int_data_sequence ,const Int_Str& int_genomic_sequence , double score_threshold , bool best_only , int min_offset , int max_offset){



	Int_Str int_data_sequence_copy = int_data_sequence;
	Int_Str int_genomic_sequence_copy = int_genomic_sequence;
	int offset_change = 0;

	if(flip_seqs){
		reverse(int_data_sequence_copy.begin(),int_data_sequence_copy.end());
		reverse(int_genomic_sequence_copy.begin(),int_genomic_sequence_copy.end());
	}

	/*if(min_offset<0){
		//Remove nucleotides that cannot be in the alignment
		//This is not true because of in/dels
		if(int_genomic_sequence.size()>(-min_offset + int_data_sequence.size())){
			short_int_genomic_sequence.erase((-min_offset + int_data_sequence.size()),string::npos);
		}
	}
	else{

	}*/
	/*if(max_offset<0){
		int_genomic_sequence_copy.erase(0,-max_offset);
		offset_change = max_offset;
	}
	else{

	}
*/
	int n_rows = int_data_sequence_copy.size();
	int n_cols = int_genomic_sequence_copy.size();



	Matrix<double> score_matrix(n_rows,n_cols);
	Matrix<int> col_memory_matrix(n_rows,n_cols);
	Matrix<int> row_memory_matrix(n_rows,n_cols);
	Matrix<int> alignment_numb_tracker(n_rows,n_cols);

	for(int i = 0 ; i!=n_rows ; ++i){
		score_matrix(i,0) = 0;
		col_memory_matrix(i,0) = 0;
		row_memory_matrix(i,0) = 0;
		for(int j = 0 ; j!=n_cols ; ++j){
			alignment_numb_tracker(i,j) = -1;
		}
	}
	for(int j = 0 ; j!=n_cols ; ++j){
		score_matrix(0,j) = 0;
		col_memory_matrix(0,j) = 0;
		row_memory_matrix(0,j) = 0;
	}

/*	int max_row_coord;
	int max_col_coord;
	int max_score = 0;*/

	vector<int> max_row_coord;
	vector<int> max_col_coord;
	vector<int> max_score;

	bool matrix_complete = false;

	int explored_row_coord = 1;
	int explored_col_coord = 1;

	bool last_column_explored = false;
	bool corner_case;
	(n_rows-1  == n_cols) ? corner_case = true : corner_case = false;




	while(!matrix_complete){


		//For efficiency the score_matrix is filled by squares at first
		//once the size of the square reaches the size of one of the sequence
		//it fills the rest

		//TODO test first if the whole genomic seq has been spanned (usually shorter than the data seq??)

		//Always start at index 1 since first column and first raw are filled with 0
		if(explored_row_coord == n_rows && !last_column_explored){
			//If all the rows have been explored
			for (int i=1 ; i!=n_rows ; ++i){
				//Explore next missing column
				sw_align_common(int_data_sequence_copy , int_genomic_sequence_copy , i , explored_col_coord-1  , score_matrix , row_memory_matrix , col_memory_matrix , alignment_numb_tracker , max_score , max_row_coord , max_col_coord);
			}

		}
		else if(explored_col_coord == n_cols){
			//If all colmuns have been explored
			for (int j=1 ; j!=n_cols ; ++j){
				//Explore next missing row
				sw_align_common(int_data_sequence_copy , int_genomic_sequence_copy , explored_row_coord-1 , j  , score_matrix , row_memory_matrix , col_memory_matrix , alignment_numb_tracker , max_score , max_row_coord , max_col_coord);
			}
			if(!last_column_explored){last_column_explored = true;}//By construction
		}
		else{
			int i=1;
			int j=1;

			while( (i!=explored_row_coord) & (j!=explored_col_coord) ){
				sw_align_common(int_data_sequence_copy , int_genomic_sequence_copy , i , explored_col_coord  , score_matrix , row_memory_matrix , col_memory_matrix , alignment_numb_tracker , max_score , max_row_coord , max_col_coord);
				++i;
				sw_align_common(int_data_sequence_copy , int_genomic_sequence_copy , explored_row_coord , j  , score_matrix , row_memory_matrix , col_memory_matrix , alignment_numb_tracker , max_score , max_row_coord , max_col_coord);
				++j;
			}
			//Fill last angle of the square
			sw_align_common(int_data_sequence_copy , int_genomic_sequence_copy , explored_row_coord , explored_col_coord  , score_matrix , row_memory_matrix , col_memory_matrix , alignment_numb_tracker , max_score , max_row_coord , max_col_coord);
		}


		if( (explored_row_coord==n_rows) & (explored_col_coord==n_cols) ){matrix_complete=true;}
		if(explored_row_coord!=n_rows){++explored_row_coord;}
		if(explored_col_coord!=n_cols){++explored_col_coord;}

	}

	//Reconstruct the alignments
	list<pair<int,Alignment_data>> seq_alignments_results;
	double max_align_score = 0;
	/*for(size_t align = 0 ; align!=max_score.size() ; align++){
		if(max_score[align]>max_align_score){
			max_align_score = max_score[align];
		}
	}*/
	for(size_t align = 0 ; align!=max_score.size() ; ++align){
		if(max_score[align]>=score_threshold){
			//if( (!best_only) | (max_score[align]==max_align_score) ){


				forward_list<int> insertions = *(new forward_list<int>); //TODO check this new
				forward_list<int> deletions = *(new forward_list<int>);
				size_t align_length = 0;

				bool end_of_alignment = false;

				int i=max_row_coord[align];
				int j=max_col_coord[align];

				size_t end_align_offset = i;

				//If sequence has been flip compute how the offset and insertion/deletion should be changed
				int flip_factor;
				int flip_offset;
				int flip_mis;
				if(flip_seqs){
					flip_factor = -1;
					flip_mis = 1;
					flip_offset = int_data_sequence.size() - int_genomic_sequence.size();
				}
				else{
					flip_factor = 1;
					flip_offset = 0;
					flip_mis = 0;
				}

			/*
				cout<<"int_data_seq: "<<int_data_sequence<<endl;
				cout<<"int_genomicseq: "<<int_genomic_sequence<<endl;
				cout<<score_matrix(max_row_coord,max_col_coord)<<endl;

					ofstream temp1(string("/home/quentin/Desktop/output_test/seq_gen_check/alignment_check/output_matrices/score_matrix.csv"));
					for(int a = 0 ; a!=n_rows ; a++){
						for(int b =0 ; b!=n_cols ; b++){
							temp1<<score_matrix(a,b)<<";";
						}
						temp1<<endl;
					}

					ofstream temp2(string("/home/quentin/Desktop/output_test/seq_gen_check/alignment_check/output_matrices/row_memory_matrix.csv"));
					for(int a = 0 ; a!=n_rows ; a++){
						for(int b =0 ; b!=n_cols ; b++){
							temp2<<row_memory_matrix(a,b)<<";";
						}
						temp2<<endl;
					}

					ofstream temp3(string("/home/quentin/Desktop/output_test/seq_gen_check/alignment_check/output_matrices/col_memory_matrix.csv"));
					for(int a = 0 ; a!=n_rows ; a++){
						for(int b =0 ; b!=n_cols ; b++){
							temp3<<col_memory_matrix(a,b)<<";";
						}
						temp3<<endl;
					}
			*/
				//TODO correct this to get the alignment until the end (not just until the best scoring nucl)

				while(!end_of_alignment){


					if( (row_memory_matrix(i,j)==0) && (col_memory_matrix(i,j)==0) ){end_of_alignment = true;}
					else if(row_memory_matrix(i,j) == 0 ){deletions.push_front(flip_factor*j+flip_mis*int_genomic_sequence_copy.size());} //TODO check the behavior of this and how to handle in-dels
					else if(col_memory_matrix(i,j) == 0 ){insertions.push_front(flip_factor*i+flip_mis*int_data_sequence_copy.size());}
					int i_temp = i;
					i-=row_memory_matrix(i_temp,j);
					j-=col_memory_matrix(i_temp,j);
					++align_length;

				}

				size_t begin_align_offset = flip_factor*i+flip_mis*int_data_sequence_copy.size();
				end_align_offset = flip_factor*end_align_offset+flip_mis*int_data_sequence_copy.size();


				//Offset is the place where the first letter of the genomic sequence aligns
				//if the alignment does not start from the beginning need to extrapolate
				int offset = flip_factor*(i-j) + flip_offset + offset_change;



				if( (offset>=min_offset) & (offset<=max_offset)){ //TODO reduce computation time by truncating the alignment from the beginning?
					//TODO change this and use incorporate_in_dels(), should probably change the list inside alignment data also to have the actual corresponding indices
					//TODO return the actual inserted/deleted sequences in the alignment data??
					Int_Str dat_seq;
					Int_Str gen_seq;
					vector<int> mismatches;
					bool neg_offset = offset<0;
					size_t n_del = distance(deletions.begin(),deletions.end());
					size_t n_ins = distance(insertions.begin(),insertions.end());
					if(neg_offset){
						gen_seq = int_genomic_sequence.substr(-offset,Int_Str::npos);
						dat_seq = int_data_sequence.substr(0,gen_seq.size()+n_ins);
					}
					else{
						dat_seq = int_data_sequence.substr(offset,Int_Str::npos);
						gen_seq = int_genomic_sequence;

					}

					if((dat_seq.size()+n_del)>(gen_seq.size()+n_ins)){
						dat_seq = dat_seq.substr(0,gen_seq.size()+n_ins-n_del);
					}
					else{
						gen_seq = gen_seq.substr(0,dat_seq.size()+n_del-n_ins);
					}

/*					cout<<"---------------------------------"<<endl;
					cout<<dat_seq.size()+n_del<<endl;
					cout<<gen_seq.size()+n_ins<<endl;*/
					size_t dat_ind=0;
					size_t gen_ind=0;

					while(dat_ind!=dat_seq.size()){

							if(neg_offset ){
								if(count(deletions.begin(),deletions.end(),gen_ind-offset)!=0){
									//The considered genomic nucleotide is deleted
									++gen_ind;
								}
								else if(count(insertions.begin(),insertions.end(),dat_ind)!=0){
									//The considered data nucleotide is an insertion
									++dat_ind;
								}
								else{
									if(not (comp_nt_int(gen_seq.at(gen_ind),dat_seq.at(dat_ind)))){
										mismatches.emplace_back(dat_ind);
									}
									++dat_ind;
									++gen_ind;
								}
							}
							else{
								if(count(deletions.begin(),deletions.end(),gen_ind)!=0){
									//The considered genomic nucleotide is deleted
									++gen_ind;
								}
								else if(count(insertions.begin(),insertions.end(),dat_ind+offset)!=0){
									//The considered data nucleotide is an insertion
									++dat_ind;
								}
								else{
									if((gen_seq.at(gen_ind)!=dat_seq.at(dat_ind))){
										mismatches.emplace_back(dat_ind+offset);
									}
									++dat_ind;
									++gen_ind;
								}
							}


					}
					if(max_score[align]>max_align_score){
						max_align_score = max_score[align];
					}
					seq_alignments_results.emplace_back(pair<int,Alignment_data>(max_score[align] , Alignment_data(offset , begin_align_offset , end_align_offset , align_length , insertions , deletions , mismatches , max_score[align])));


				}

		}
	}
	if(best_only & (seq_alignments_results.size()>1)){
		for(list<pair<int,Alignment_data>>::iterator align = seq_alignments_results.begin() ; align != seq_alignments_results.end() ; ++align ){
			if((*align).first < max_align_score){
				align = seq_alignments_results.erase(align);
				--align;
			}
		}
	}



	return seq_alignments_results;
}

//Compare alignments (sort by score)
bool align_compare(Alignment_data align1 , Alignment_data align2 ){
	return align1.score > align2.score;
}

/**
 * A dumb function to read CSV anchor gene indices
 */
unordered_map<string,size_t> read_gene_anchors_csv(string filename , string sep){
	ifstream infile(filename);
		if(!infile){
			throw runtime_error("File not found: "+filename + " in read_gene_anchors_csv()");
		}

		string temp_str;
		unordered_map<string,size_t> anchors_map;

		getline(infile,temp_str); //Ignore first header line

		while (getline(infile,temp_str)){
			int coma_index = temp_str.find(sep);
			anchors_map.emplace(temp_str.substr(0,coma_index),stoi(temp_str.substr(coma_index+1,string::npos)));
		}

		return anchors_map;
}


