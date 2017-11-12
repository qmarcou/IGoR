import pandas as pd
from ..utils.utils import get_str_asarray

def compute_mutation_frequency(x,threshold): 
	"""Computes mutation frequency for positions with coverage larger than threshold."""
    fraction_list = []
    for cov,err in zip(x.coverage,x.errors):
        if(cov>threshold):
            fraction_list.append(err/cov)
        else:
            fraction_list.append(np.NaN)
    return fraction_list


#Read gene coverage and errors data
def read_coverage_and_errors_file(filename,get_diag_of_N_dim=None):
	"""Reads coverage and error counter file as pandas.DataFrame

	For Ndimensional joint coverage and errors the diagonal(equivalent to 1D) can be extracted by passing the number of dimensions as an argument.
	"""
    raw_read = pd.read_csv(filename,sep=';')
    
    #Convert coverage strings into arrays
    tmp_cov = raw_read.apply(lambda x : get_str_asarray(x.coverage),axis=1)
    #Convert errors strings into arrays
    tmp_err = raw_read.apply(lambda x : get_str_asarray(x.errors),axis=1)
    
    if get_diag_of_N_dim != None:
        for i in range(0,len(tmp_cov)):
            single_dim_len = int(float(len(tmp_cov.iloc[i]))**(1.0/float(get_diag_of_N_dim)))
            reshape_tuple = tuple([single_dim_len]*get_diag_of_N_dim)
            #print(reshape_tuple)
            #print(len(tmp_cov))

            #Get the diagonal
            tmp_cov.iloc[i] = list(np.asarray(tmp_cov.iloc[i]).reshape(reshape_tuple).diagonal())
            tmp_err.iloc[i] = list(np.asarray(tmp_err.iloc[i]).reshape(reshape_tuple).diagonal())
    
    raw_read.coverage = tmp_cov
    raw_read.errors = tmp_err
    
    return raw_read
