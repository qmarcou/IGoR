import pandas

def extract_best_aligns(aligns_df):
	mask=aligns_df.groupby('seq_index').agg({'score':'idxmax'}) #get /!\ FIRST /!\ index of max align score for each sequence
	aligns_df_best= aligns_df.loc[mask['score']].reset_index(drop=True)
	return aligns_df_best

def get_misinsdel_asarray(misinsdel_str):
    misinsdel_array = []
    if(len(misinsdel_str)>2):
        misinsdel_str = misinsdel_str[1:-1]
        next_comma_index = misinsdel_str.find(',')
        comma_index = -1
        while next_comma_index!=-1:
             misinsdel_array.append(int(misinsdel_str[comma_index+1:next_comma_index]))
             comma_index = next_comma_index
             next_comma_index = misinsdel_str.find(',',next_comma_index+1)
        misinsdel_array.append(int(misinsdel_str[comma_index+1:]))     
    return misinsdel_array


def read_alignments(filename):
	aligns = pandas.read_csv(filename,delimiter=';')
	#Convert the string of insertions into an array of integers
	tmp = aligns.apply(lambda x: get_misinsdel_asarray(x.insertions),axis=1);
	aligns.insertions = tmp
	#Convert the string of deletions into an array of integers
	tmp = aligns.apply(lambda x: get_misinsdel_asarray(x.deletions),axis=1);
	aligns.deletions = tmp
	#Convert the string of mismatches into an array of integers
	tmp = aligns.apply(lambda x: get_misinsdel_asarray(x.mismatches),axis=1);
	aligns.mismatches = tmp
	
	return aligns

def read_best_alignments(filename):
	#Not efficient just faster to code

	aligns = pandas.read_csv(filename,delimiter=';')
	aligns = extract_best_aligns(aligns)

	#Convert the string of insertions into an array of integers
	tmp = aligns.insertions.apply(get_misinsdel_asarray);
	aligns.insertions = tmp
	#Convert the string of deletions into an array of integers
	tmp = aligns.deletions.apply(get_misinsdel_asarray);
	aligns.deletions = tmp
	#Convert the string of mismatches into an array of integers
	tmp = aligns.mismatches.apply(get_misinsdel_asarray);
	aligns.mismatches = tmp
	
	return aligns

#Import genomic template sequences
import re as regex
def read_FASTA_strings(filename):
	seq = regex.compile('>')
	line = regex.compile('\n')
	with open(filename) as file:
		tmp = seq.split(file.read())
		final = {}
		for i in range(1,len(tmp)):
			split_seq = line.split(tmp[i])
			final[split_seq[0]]=split_seq[1].upper()
		return final

	
def has_mismatches(x):
    tmp = asarray(x.mismatches)
    return any(multiply(tmp>= x["5_p_align_offset"],tmp<=x["3_p_align_offset"]))	
