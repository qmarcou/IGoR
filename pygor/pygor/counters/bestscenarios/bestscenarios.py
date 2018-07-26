#      Author: Quentin Marcou
#
#  This source code is distributed as part of the IGoR software.
#  IGoR (Inference and Generation of Repertoires) is a versatile software to
#  analyze and model immune receptors generation, selection, mutation and all
#  other processes.
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

# Declare functions translating scenarios indices into "real values"

import copy

import pandas

from ...models.genmodel import GenModel
from ...utils import utils


def read_bestscenarios_values(scenarios_file, model_file):
    best_scenarios_df = read_bestscenarios_indices(scenarios_file)
    genmodel = GenModel(model_file)
    return scenarios_indices2values(best_scenarios_df, genmodel)


def read_bestscenarios_indices(filename):
    return pandas.read_csv(filename, sep=';')


def scenarios_indices2values(best_scenarios, input_genmodel,
                             drop_alleles=False, name2nickname=True):
    """ This function uses a supplied generative model to transform realization
    indices into realization values.

    """
    best_scenarios_real = copy.deepcopy(best_scenarios)
    for event in input_genmodel.events:
        # Enable truncated scenarios_dataframe to be processed
        if list(best_scenarios_real.keys()).count(event.name) > 0:
            real_vect = event.get_realization_vector()

            # Now get realizations as integers or list of integers
            if event.event_type == "DinucMarkov":
                tmp_real_indices = best_scenarios[event.name].apply(
                    lambda x: utils.get_str_asarray(x, dtype=int,
                                                    boundaries_char=["(", ")"],
                                                    sep=','))

            else:
                tmp_real_indices = best_scenarios[event.name].apply(
                    get_real_str_as_int)
            best_scenarios_real[event.name] = tmp_real_indices.apply(
                lambda x: real_vect[x])

            # Return mismatches as a list of positions
            best_scenarios_real['Mismatches'] = best_scenarios[
                'Mismatches'].apply(lambda x: utils.get_str_asarray(
                    x, dtype=int, boundaries_char=["(", ")"], sep=','))

            if drop_alleles and (event.event_type == "GeneChoice"):
                best_scenarios_real[event.name] = best_scenarios_real[
                    event.name].apply(drop_allele_info)

            if name2nickname:
                best_scenarios_real.rename(
                    columns={event.name: event.nickname}, inplace=True)

    return best_scenarios_real


def drop_allele_info(x):
    if type(x) == str:
        return x[0:x.find('*')]
    else:
        return None


def get_real_str_as_int(real_str):
    if type(real_str) == str:
        return int(real_str[1:-1])
    else:
        return None
