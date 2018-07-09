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

import copy

import numpy
from pylab import find


def compute_cross_entropy_subparts(event_name, model1, model2,
                                   debug_output=False):
    """Returns the cross entropy -sum(p_model1*log(p_model2)) for the specified
    recombination event.

    """
    # Compute -sum(p_model1*log(p_model2)) for one element of the model

    # TODO make some background check on model_parms

    for event in model1.events:
        if event.name == event_name:
            event_cluster = list()
            event_cluster.append(event)

            # Get all events related (even remotely) to this considered event
            explored = False
            parents_names = list()
            parents_names.append(event.name)
            while not explored:
                next_parents = list()
                for parent in parents_names:
                    if parent not in model1.edges:
                        continue
                    parent_parents = model1.edges[parent].parents
                    for p in parent_parents:
                        next_parents.append(p)
                        if event_cluster.count(model1.get_event(p)) == 0:
                            event_cluster.append(model1.get_event(p))
                parents_names = next_parents
                if len(parents_names) == 0:
                    explored = True

            dimension_list = []
            dimension_names_list = []
            for ev in event_cluster:
                dimension_list.append(
                    model1.marginals[0][ev.nickname].shape[-1])
                dimension_names_list.append(ev.nickname)

            if event.event_type != "DinucMarkov":

                final_array = numpy.ones(tuple(dimension_list))
                # print(final_array.shape)

                # Reshape all marginals in one big redundant array
                reshaped_marginals = list()

                # Takes care of model 1(non log part)
                for ev in event_cluster:
                    if ev.name not in model1.edges:
                        event_parents = []
                    else:
                        event_parents = model1.edges[ev.name].parents
                    # print(event_parents)
                    dimension_list = []
                    for e in event_cluster:
                        if ((event_parents.count(e.name) > 0)
                                or (ev.name == e.name)):
                            dimension_list.append(
                                model1.marginals[0][e.nickname].shape[-1])  # same as getting the number of realizations
                        else:
                            dimension_list.append(1)

                        # print(ev)
                        # print(dimension_list)
                        # print(shape(model1.marginals[ev.nickname]))
                        # print("Dimension list")
                        # print(dimension_list)

                    swapped_marginal = copy.deepcopy(
                        model1.marginals[0][ev.nickname])
                    swapped_dimension_names = copy.deepcopy(
                        model1.marginals[1][ev.nickname])

                    is_sorted = False
                    while not is_sorted:
                        for i in range(0, len(swapped_marginal.shape)):
                            if i == len(swapped_marginal.shape) - 1:
                                is_sorted = True
                                break
                            elif (find(numpy.asarray(dimension_names_list) ==
                                       swapped_dimension_names[i]) > find(
                                    numpy.asarray(dimension_names_list) ==
                                    swapped_dimension_names[i + 1])):
                                swapped_marginal = swapped_marginal.swapaxes(
                                    i, i + 1)
                                # Artisanal swap
                                tmp = swapped_dimension_names[i + 1]
                                swapped_dimension_names[i + 1] = \
                                    swapped_dimension_names[i]
                                swapped_dimension_names[i] = tmp
                                break

                    reshaped_marginals.append(
                        swapped_marginal.reshape(tuple(dimension_list)))

                # Takes care of model 2(the log part)
                # event_parents = model1.edges[event_name].parents
                if event_name not in model1.edges:
                    event_parents = []
                else:
                    event_parents = model1.edges[event_name].parents
                dimension_list = []
                for e in event_cluster:
                    if ((event_parents.count(e.name) > 0)
                            or (event_name == e.name)):
                        dimension_list.append(
                            model1.marginals[0][e.nickname].shape[-1])  # same as getting the number of realizations
                    else:
                        dimension_list.append(1)

                    swapped_marginal = copy.deepcopy(
                        model2.marginals[0][event.nickname])
                    swapped_dimension_names = copy.deepcopy(
                        model2.marginals[1][event.nickname])

                    is_sorted = False
                    while not is_sorted:
                        # print(swapped_marginal.shape)
                        for i in range(0, len(swapped_marginal.shape)):
                            if i == len(swapped_marginal.shape) - 1:
                                is_sorted = True
                                break
                            elif (find(numpy.asarray(dimension_names_list) ==
                                       swapped_dimension_names[i]) > find(
                                    numpy.asarray(dimension_names_list) ==
                                    swapped_dimension_names[i + 1])):
                                swapped_marginal = swapped_marginal.swapaxes(
                                    i, i + 1)
                                # Artisanal swap
                                tmp = swapped_dimension_names[i + 1]
                                swapped_dimension_names[i + 1] = \
                                    swapped_dimension_names[i]
                                swapped_dimension_names[i] = tmp
                                break

                reshaped_event_marginals = swapped_marginal.reshape(
                    tuple(dimension_list))

                for reshaped_marg in reshaped_marginals:
                    final_array *= reshaped_marg
                # print(final_array)
                # print(log(reshaped_event_marginals))
                # Add pseudo-counts for 0 entries
                pseudo_count = reshaped_event_marginals[
                                   reshaped_event_marginals != 0]\
                                   .min() * 10 ** (-5)
                final_array *= numpy.log2(
                    reshaped_event_marginals + pseudo_count)
                if debug_output:
                    return -final_array
                else:
                    return -final_array.sum()

            if event.event_type == "DinucMarkov":
                # If VDJ needs to compute separate entropy contributions for VD and DJ
                junction_list = []
                if event.seq_type == "VDJ_genes":
                    junction_list = ["VD_genes", "DJ_gene"]
                else:
                    junction_list = [event.seq_type]
                # First need to find the corresponding insertion event

                dinuc_cross_entropy = []
                for junction in junction_list:
                    # First need to find the corresponding insertion event
                    insertion_event = None
                    junction_found = False
                    for ev in model1.events:
                        if ev.event_type == "Insertion" \
                                and ev.seq_type == junction:
                            insertion_event = ev
                            junction_found = True
                            break
                    if not junction_found:
                        print("INSERTION JUNCTION NOT FOUND: " + junction)
                    dinuc_cross_entropy.append((junction,
                                                compute_markov_steady_entropy(
                                                    event, model1, model2,
                                                    event_cluster,
                                                    dimension_list,
                                                    dimension_names_list,
                                                    insertion_event)))
                if len(junction_list) == 1:
                    return dinuc_cross_entropy[0][1]  # Only return a double
                else:
                    return dinuc_cross_entropy


def compute_markov_steady_entropy(dinuc_event, model1, model2, event_cluster,
                                  dimension_list, dimension_names_list,
                                  insertion_event):
    """Compute the Dinucleotide markov model cross entropy. The entropy is
    approximated assuming the markov chain is at steady state.

    """
    final_array = numpy.ones(tuple(dimension_list))
    # print(final_array.shape)

    # Reshape all marginals in one big redundant array
    reshaped_marginals = list()

    # Reverse event_cluster ? (to make sure dinuc marginals is the last dimension), or just swap it afterwards

    for ev in event_cluster:
        if ev.name not in model1.edges:
            event_parents = []
        else:
            event_parents = model1.edges[ev.name].parents
        # print(event_parents)
        dimension_list = []
        tot_dimension = 1
        for e in event_cluster:
            if event_parents.count(e.name) > 0:  # Dimension for the dinuc is 1 since we'll store vectors and matrices as dtype
                dimension_list.append(model1.marginals[0][e.nickname]
                                      .shape[-1])  # same as getting the number of realizations
                tot_dimension *= model1.marginals[0][e.nickname].shape[-1]
            else:
                dimension_list.append(1)
        # print("Tot dimension:" + str(tot_dimension))

        # Swap marginals in the correct order for model1 contributions
        swapped_marginal = copy.deepcopy(model1.marginals[0][ev.nickname])
        swapped_dimension_names = copy.deepcopy(
            model1.marginals[1][ev.nickname])

        is_sorted = False
        while not is_sorted:
            for i in range(0, len(swapped_marginal.shape)):
                if i == len(swapped_marginal.shape) - 1:
                    is_sorted = True
                    break
                elif (find(numpy.asarray(dimension_names_list) ==
                           swapped_dimension_names[i]) > find(
                        numpy.asarray(dimension_names_list) ==
                        swapped_dimension_names[i + 1])):
                    swapped_marginal = swapped_marginal.swapaxes(i, i + 1)
                    # Artisanal swap
                    tmp = swapped_dimension_names[i + 1]
                    swapped_dimension_names[i + 1] = swapped_dimension_names[i]
                    swapped_dimension_names[i] = tmp
                    break

        # DINUC SHOULD PROBABLY NOT BE IN THE EVENT CLUSTER

        if ev != dinuc_event:
            reshaped_marginals.append(
                swapped_marginal.reshape(tuple(dimension_list)))
        else:
            flattened_marginals_1 = swapped_marginal
    # print(reshaped_marginals)
    # print(flattened_marginals_1)

    # Swap model 2 marginals
    swapped_marginal = copy.deepcopy(model2.marginals[0][dinuc_event.nickname])
    swapped_dimension_names = copy.deepcopy(
        model2.marginals[1][dinuc_event.nickname])

    is_sorted = False
    while not is_sorted:
        # print(swapped_marginal.shape)
        for i in range(0, len(swapped_marginal.shape)):
            if i == len(swapped_marginal.shape) - 1:
                is_sorted = True
                break
            elif (find(numpy.asarray(dimension_names_list) ==
                       swapped_dimension_names[i]) > find(
                    numpy.asarray(dimension_names_list) ==
                    swapped_dimension_names[i + 1])):
                swapped_marginal = swapped_marginal.swapaxes(i, i + 1)
                # Artisanal swap
                tmp = swapped_dimension_names[i + 1]
                swapped_dimension_names[i + 1] = swapped_dimension_names[i]
                swapped_dimension_names[i] = tmp
                break

    # Get dinuc marginals
    # flattened_marginals_1 = copy.deepcopy(model1.marginals[0][dinuc_event.nickname])
    flattened_marginals_1 = flattened_marginals_1.reshape(
        (tot_dimension, 4, 4))

    mono_nt_vect_1 = numpy.ndarray((tot_dimension), dtype=numpy.matrix)
    conditional_nt_mat_1 = numpy.ndarray((tot_dimension), dtype=numpy.matrix)

    flattened_marginals_2 = swapped_marginal
    flattened_marginals_2 = flattened_marginals_2.reshape(
        (tot_dimension, 4, 4))

    mono_nt_vect_2 = numpy.ndarray(tot_dimension, dtype=numpy.matrix)
    conditional_nt_mat_2 = numpy.ndarray(tot_dimension, dtype=numpy.matrix)

    mono_nt_cross_entropy = numpy.ndarray(tot_dimension, dtype=float)  # entropy of steady state distribution
    cross_entropy_production_rate = numpy.ndarray(tot_dimension, dtype=float)  # entropy of steady state distribution

    for i in range(0, tot_dimension):
        transition_matrix_1 = numpy.mat(flattened_marginals_1[i, :, :])
        transition_matrix_2 = numpy.mat(flattened_marginals_2[i, :, :])
        for [transition_matrix, mono_nt_vect, conditional_nt_mat] in [
            (transition_matrix_1, mono_nt_vect_1, conditional_nt_mat_1),
            (transition_matrix_2, mono_nt_vect_2, conditional_nt_mat_2)]:
            # print(transition_matrix)
            [eigenValue, eigenVect] = numpy.linalg.eig(
                transition_matrix.transpose())

            # Get the unit eigenvalue
            unit_eig = find(abs(eigenValue - 1) < .00001)
            if len(unit_eig) > 1:
                print("Several unit eigenvalues for model1 transition matrix")

            # Get the corresponding eigenvector and make sure it is normalized (numpy returns a unit vector)
            stat_dist = eigenVect[:, unit_eig].transpose()
            if stat_dist.any() < 0:
                if stat_dist.all() <= 0:
                    stat_dist *= (-1)
                else:
                    print("POSITIVE AND NEGATIVE ENTRY IN LEFT EIGENVECTOR")
            stat_dist = numpy.divide(stat_dist, sum(
                stat_dist))  # although should already be a unit vector

            # print("Stationnary distribution:" + str(stat_dist))

            mono_nt_vect[i] = stat_dist
            conditional_nt_mat[i] = numpy.diag(
                numpy.asarray(stat_dist).reshape(4))
        # print(numpy.asarray(stat_dist).reshape(4))
        # print(numpy.diag(numpy.asarray(stat_dist).reshape((4))))

        mono_nt_cross_entropy[i] = numpy.nansum(
            numpy.multiply(mono_nt_vect_1[i], numpy.log2(
                mono_nt_vect_2[i])))  # SHOULD PROBABLY ADD PSEUDOCOUNTS
        cross_entropy_production_rate[i] = numpy.nansum(numpy.nansum(
            numpy.multiply(conditional_nt_mat_1[i] * transition_matrix_1,
                           numpy.log2(transition_matrix_2))))

    # print("Mono nt entropy:"+str(mono_nt_cross_entropy[i]))
    # print("Cross entropy production: "+str(cross_entropy_production_rate[i]))

    # Check if the dinucleotide markov model depends on the insertion event
    insertion_dependent = False
    if dinuc_event.name in model1.edges and model1.edges.parents.count(
            insertion_event.name):
        insertion_dependent = True

    cross_entropy = numpy.ndarray((tot_dimension), dtype=float)
    if insertion_dependent:
        print("INSERTION DEPENDENT DINUC MARKOV MODEL NOT HANDLED SO FAR")
    else:
        insertion_marginals = model1.marginals[0][insertion_event.nickname]
        for i in range(0, tot_dimension):
            insertion_weighed_entropy = 0
            tmp = 0
            for ins_real in insertion_event.realizations:
                L_ins_entropy = 0

                if ins_real.value > 0:
                    L_ins_entropy += mono_nt_cross_entropy[i]  # first nucleotide drawn from steady state distribution
                    L_ins_entropy += cross_entropy_production_rate * (
                                ins_real.value - 1)  # add cross entropy production for every markov step
                L_ins_entropy *= insertion_marginals[ins_real.index]
                tmp += insertion_marginals[ins_real.index] * 2 * ins_real.value
                insertion_weighed_entropy += L_ins_entropy  # MISTAKE DOES NOT TAKE INTO ACCOUNT DEPENDENCIES OF THE INSERTION DISTRIBUTION
            cross_entropy[i] = insertion_weighed_entropy
    # print("Insertion weighed cross entropy:"+str(cross_entropy))
    # print("entropy_bound: "+str(tmp))

    final_array = cross_entropy.reshape(tuple(dimension_list))
    # Now weigh the entropy by all dependences realization
    for reshaped_marg in reshaped_marginals:
        final_array *= reshaped_marg

    return -final_array.sum()


def compute_event_entropy(event_name, model, debug_output=False):
    """Compute the provided recombination event entropy."""
    return compute_cross_entropy_subparts(event_name, model, model,
                                          debug_output)


def compute_event_DKL(event_name, model1, model2, debug_output=False):
    """Compute the Kullback Leibler divergence for the specified event and
    provided models. Returns  D_KL(model1||model2) for the considered event.

    """
    return compute_cross_entropy_subparts(event_name, model1, model2,
                                          debug_output) - \
        compute_event_entropy(event_name, model1, debug_output)


def compute_models_DKL(model1, model2):
    """Compute the Kullback Leibler divergence two recombination models.
    Returns D_KL(model1||model2).

    """
    total_DKL = 0
    for event in model1.events:
        total_DKL += compute_event_DKL(event.name, model1, model2)
    return total_DKL


def compute_model_entropy(model):
    """Compute the recombination entropy of the supplied GenModel."""
    total_entropy = 0
    for event in model.events:
        total_entropy += compute_event_entropy(event.name, model)
    return total_entropy


def compute_type_model_DKL(type_list, model1, model2):
    """Compute Kullback Leibler divergence associated with events of supplied
    types.

    """
    total_DKL = 0
    for event in model1.events:
        if (numpy.asarray(type_list) == event.event_type).any():
            total_DKL += compute_event_DKL(event.name, model1, model2)
    return total_DKL
