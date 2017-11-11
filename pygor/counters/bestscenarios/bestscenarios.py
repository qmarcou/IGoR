#Declare functions translating scenarios indices into "real values"
from ..models import genmodel 
from ..utils import utils
import pandas as pd
import copy

def read_bestscenarios_indices(filename):
	return pd.read_csv(filename,sep=';')

def read_bestscenarios_values(scenarios_file,model_file):
	best_scenarios_df = read_bestscenarios_values(scenarios_file)
	genmodel = genmodel.GenModel(model_file)
	return scenarios_indices2values(best_scenarios_df,genmodel)

def scenarios_indices2values(best_scenarios,input_genmodel,drop_alleles=False,name2nickname=True):
#""" This function uses a supplied generative model to transform realization indices into realization values """ 
	best_scenarios_real = copy.deepcopy(best_scenarios)
	for event in input_genmodel.events:
		# Enable truncated scenarios_dataframe to be processed
		if list(best_scenarios_real.keys()).count(event.name)>0:
			real_vect = event.get_realization_vector()
	
			# Now get realizations as integers or list of integers
			if event.event_type == "DinucMarkov":
				tmp_real_indices = best_scenarios[event.name].apply(lambda x:get_str_asarray(x,dtype=int,boundaries_char = ["(",")"],sep = ','))

			else:
				tmp_real_indices = best_scenarios[event.name].apply(get_real_str_as_int)
			best_scenarios_real[event.name] = tmp_real_indices.apply(lambda x: real_vect[x])

			#Return mismatches as a list of positions
			best_scenarios_real['Mismatches'] = best_scenarios['Mismatches'].apply(lambda x:get_str_asarray(x,dtype=int,boundaries_char = ["(",")"],sep = ','))


			if drop_alleles and (event.event_type == "GeneChoice"):
				best_scenarios_real[event.name] = best_scenarios_real[event.name].apply(drop_allele_info)
	
			if name2nickname:
				best_scenarios_real.rename(columns={event.name:event.nickname},inplace=True)			

	return best_scenarios_real


def drop_allele_info(x):
    if type(x)==str:
        return x[0:x.find('*')]
    else:
        return None

def get_real_str_as_int(real_str):
    if type(real_str)==str:
        return int(real_str[1:-1])
    else:
        return None


