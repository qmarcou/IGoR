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
