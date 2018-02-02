 #      Author: Quentin Marcou
 #      
 #  This source code is distributed as part of the IGoR software.
 #  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 #  generation, selection, mutation and all other processes.
 #   Copyright (C) 2017  Quentin Marcou
 #
 #   This program is free software: you can redistribute it and/or modify
 #   it under the terms of the GNU General Public License as published by
 #   the Free Software Foundation, either version 3 of the License, or
 #   (at your option) any later version.
 #
 #   This program is distributed in the hope that it will be useful,
 #   but WITHOUT ANY WARRANTY; without even the implied warranty of
 #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #   GNU General Public License for more details.

 #   You should have received a copy of the GNU General Public License
 #   along with this program.  If not, see <https://www.gnu.org/licenses/>.

def nt2int(seq):
    int_seq = []
    for nt in seq:
        if nt=='A':
            int_seq.append(0)
        elif nt=='C':
            int_seq.append(1)
        elif nt=='G':
            int_seq.append(2)
        else:
            int_seq.append(3)
    return int_seq

def int2nt(seq):
    nt_seq = []
    for int_nt in seq:
        if int_nt==0:
            nt_seq.append('A')
        elif int_nt==1:
            nt_seq.append('C')
        elif int_nt==2:
            nt_seq.append('G')
        else:
            nt_seq.append('T')
    return nt_seq

def get_seq_rev_comp(seq):
    if not type(seq)==str:
        print("Must provide nucleotide sequence for get seq rev comp function")
    rev_comp_seq = ""
    seq = seq.upper()
    for nt in seq:
        if nt=="A":
            rev_comp_seq = "T" + rev_comp_seq
        elif nt=="T":
            rev_comp_seq = "A" + rev_comp_seq
        elif nt=="C":
            rev_comp_seq = "G" + rev_comp_seq
        elif nt=="G":
            rev_comp_seq = "C" + rev_comp_seq
        else:
            print("Nucleotide " + nt + " not supported in get seq rev comp function")
    return rev_comp_seq

def get_str_asarray(full_str,dtype=float,boundaries_char = ["(",")"],sep = ','):
#""" Returns subfields of string seperated by *sep* , remove boundary chars and convert to the provided type """
    full_array = []
    if(len(full_str)>(len(boundaries_char[0])+len(boundaries_char[1]))):
        
        #Check if the beginning and end of the string indeed match the boundaries characters 
        if full_str.find(boundaries_char[0])!=0:
            print("Beginning char cannot be found")
        elif full_str.find(boundaries_char[1])!=(len(full_str)-len(boundaries_char[1])):
            print("Ending char cannot be found")
            
        full_str = full_str[len(boundaries_char[0]):-len(boundaries_char[1])]
        next_comma_index = full_str.find(sep)
        comma_index = -1
        while next_comma_index!=-1:
             full_array.append(dtype(full_str[comma_index+1:next_comma_index]))
             comma_index = next_comma_index
             next_comma_index = full_str.find(',',next_comma_index+1)
        full_array.append(dtype(full_str[comma_index+1:]))     
    return full_array
