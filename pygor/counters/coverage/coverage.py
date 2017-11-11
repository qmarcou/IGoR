import pandas as pd
from ..utils import utils

def compute_mutation_frequency(x,threshold=2000.0): #=800
    fraction_list = []
    for cov,err in zip(x.coverage,x.errors):
        if(cov>threshold):
            fraction_list.append(err/cov)
        else:
            fraction_list.append(np.NaN)
    return fraction_list



def get_str_asarray(full_str,dtype=float,boundaries_char = ["(",")"],sep = ','):
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
        full_array.append(float(full_str[comma_index+1:]))     
    return full_array

#Read gene coverage and errors data
def read_coverage_and_errors_file(filename,get_diag_of_N_dim=None):
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
