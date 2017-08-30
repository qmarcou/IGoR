#Declare functions translating scenarios indices into "real values"
from ..models import genmodel 
import pandas as pd

def read_bestscenarios_indices(filename):
	return best_scenarios_df = pd.read_csv(filename,sep=';')

def read_bestscenarios_values(scenarios_file,model_file):
	best_scenarios_df = read_bestscenarios_values(scenarios_file)
	genmodel = genmodel.GenModel(model_file)
	return scenarios_indices2values(best_scenarios_df,genmodel)

def scenarios_indices2values(best_scenarios,genmodel,old_columns_list,new_columns_names):
	""" Write some documentation """ 
    best_scenarios_real = pd.DataFrame()
    
    for column_name,new_column_name in zip(old_columns_list,new_columns_names):
        #Get the correct event
        event = genmodel.get_event(column_name)
        real_vect = event.get_realization_vector()

        #Deal with best scenarios
        #Now get realizations as integers
        tmp_real_indices = best_scenarios[column_name].apply(get_real_str_as_int)
        best_scenarios_real[new_column_name] = tmp_real_indices.apply(lambda x: real_vect[x])

    best_scenarios_real['seq_index'] = best_scenarios['seq_index']
    if any(best_scenarios.columns=='scenario_proba_cond_seq'):
        best_scenarios_real['scenario_proba_cond_seq'] = best_scenarios['scenario_proba_cond_seq']
    best_scenarios_real['scenario_rank'] = best_scenarios['scenario_rank']
    best_scenarios_real = best_scenarios_real[best_scenarios_real.seq_index.notnull()]


    for gene_column in ['V_gene','D_gene','J_gene']:
        best_scenarios_real[gene_column] = best_scenarios_real[gene_column].apply(drop_allele_info)

    #for best_scenarios_real,best_scenarios in zip([best_scenarios_viterbi_real,best_scenarios_learned_real,best_scenarios_true_real],[best_scenarios_viterbi,best_scenarios_learned,best_scenarios_true]):
        #best_scenarios_real = best_scenarios_real.reindex_like(sorted(best_scenarios_real.columns),axis=1,copy=False)
    return best_scenarios_real

#Declare a function taking the true realizations and translating them into "real values"

def true_realizations_indices2values(true_realizations,genmodel,old_columns_list,new_columns_names):
    true_realizations_real = pd.DataFrame()

    for column_name,new_column_name in zip(old_columns_list,new_columns_names):
        #Get the correct event
        event = genmodel.get_event(column_name)
        real_vect = event.get_realization_vector()

        #Deal with true realizations
        #Now get realizations as integers
        tmp_real_indices = true_realizations[column_name].apply(get_real_str_as_int)
        true_realizations_real[new_column_name] = tmp_real_indices.apply(lambda x: real_vect[x])


    for gene_column in ['V_gene','D_gene','J_gene']:
        true_realizations_real[gene_column] = true_realizations_real[gene_column].apply(drop_allele_info)

    true_realizations_real = true_realizations_real.reindex_axis(sorted(true_realizations_real.columns), axis=1, copy=False);
    return true_realizations_real

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
