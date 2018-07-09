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


class GenModel:
    """Generative model class containing model topology and realization
    probabilities. Similar to IGoR's C++ GenModel class.

    """
    def __init__(self):
        self.error_rate = numpy.NaN
        self.events = list()
        self.edges = {}  # Empty dictionary
        self.marginals = {}

    def __init__(self):
        self.error_rate = numpy.NaN
        self.events = list()
        self.edges = {}  # Empty dictionary
        self.marginals = {}

    def __init__(self, model_parms_file):
        self.error_rate = numpy.NaN
        self.events = list()
        self.edges = {}  # Empty dictionary
        self.marginals = {}
        self.read_model_parms(model_parms_file)
        self.marginals = {}

    def __init__(self, model_parms_file, marginals_file):
        self.error_rate = numpy.NaN
        self.events = list()
        self.edges = {}  # Empty dictionary
        self.marginals = {}
        self.read_model_parms(model_parms_file)
        self.marginals = read_marginals_txt(marginals_file)

    def read_model_parms(self, filename):
        """Reads a model graph structure from a model params file.
        Note that for now this method does not read the error rate information.

        """
        with open(filename, 'r') as file:
            # dictionary containing recarrays?
            value_str = []
            index = []
            real_name = []
            current_part = ""

            line = file.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip(
                '\r')  # Remove carriage return character (if needed)

            if strip_line == "@Event_list":
                line = file.readline()
                strip_line = line.rstrip('\n')  # Remove end of line character
                strip_line = strip_line.rstrip(
                    '\r')  # Remove carriage return character (if needed)

                while strip_line[0] == '#':
                    semicolon_index = strip_line.find(';')
                    event_type = strip_line[1:semicolon_index]

                    next_semicolon_index = strip_line.find(';',
                                                           semicolon_index + 1)
                    seq_type = strip_line[
                               semicolon_index + 1:next_semicolon_index]

                    semicolon_index = next_semicolon_index
                    next_semicolon_index = strip_line.find(';',
                                                           semicolon_index + 1)
                    seq_side = strip_line[
                               semicolon_index + 1:next_semicolon_index]

                    semicolon_index = next_semicolon_index
                    next_semicolon_index = strip_line.find(';',
                                                           semicolon_index + 1)
                    priority = strip_line[
                               semicolon_index + 1:next_semicolon_index]

                    semicolon_index = next_semicolon_index
                    next_semicolon_index = strip_line.find(';',
                                                           semicolon_index + 1)
                    nickname = strip_line[semicolon_index + 1:]

                    event = Rec_Event(event_type, seq_type, seq_side, priority,
                                      nickname)

                    line = file.readline()
                    strip_line = line.rstrip('\n')  # Remove end of line character
                    strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
                    while strip_line[0] == '%':
                        name = ""
                        semicolon_index = 0
                        next_semicolon_index = strip_line.find(';')
                        if event_type == "GeneChoice":
                            name = strip_line[
                                   semicolon_index + 1:next_semicolon_index]
                            semicolon_index = next_semicolon_index
                            next_semicolon_index = \
                                strip_line.find(';', semicolon_index + 1)
                            value = strip_line[
                                    semicolon_index + 1:next_semicolon_index]
                        elif event_type == "DinucMarkov":
                            value = strip_line[
                                    semicolon_index + 1:next_semicolon_index]
                        else:
                            value = int(strip_line[
                                        semicolon_index + 1:next_semicolon_index])

                        semicolon_index = next_semicolon_index
                        next_semicolon_index = strip_line\
                            .find(';', semicolon_index + 1)

                        index = int(strip_line[semicolon_index + 1:])

                        realization = Event_realization(name, value, index)

                        event.add_realization(realization)

                        line = file.readline()
                        strip_line = line.rstrip('\n')  # Remove end of line character
                        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
                    self.events.append(event)

            if strip_line == "@Edges":

                line = file.readline()
                strip_line = line.rstrip('\n')  # Remove end of line character
                strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
                while strip_line[0] == '%':
                    semicolon_index = strip_line.find(';')
                    parent = strip_line[1:semicolon_index]
                    child = strip_line[semicolon_index + 1:]

                    if not (parent in self.edges):
                        self.edges[parent] = Adjacency_List()

                    if not (child in self.edges):
                        self.edges[child] = Adjacency_List()

                    self.edges[parent].children.append(child)
                    self.edges[child].parents.append(parent)

                    line = file.readline()
                    strip_line = line.rstrip('\n')  # Remove end of line character
                    strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

    def get_event(self, event_name, by_nickname=False):
        """Returns the RecEvent with corresponding name or nickname."""
        if by_nickname:
            for ev in self.events:
                if ev.nickname == event_name:
                    return ev
            raise Exception(
                'RecEvent with nickname \"' + event_name + "\" not found.")
        else:
            for ev in self.events:
                if ev.name == event_name:
                    return ev
            raise Exception(
                'RecEvent with name \"' + event_name + "\" not found.")


class Rec_Event:
    """Recombination event class containing event's name, type, realizations,
    etc... Similar to IGoR's C++ RecEvent class.

    """
    def __init__(self, event_type, seq_type, seq_side, priority):
        self.event_type = event_type
        self.seq_type = seq_type
        self.seq_side = seq_side
        self.priority = priority
        self.realizations = list()
        self.name = ""
        self.nickname = ""
        self.update_name()
        return

    def __init__(self, event_type, seq_type, seq_side, priority, nickname):
        self.event_type = event_type
        self.seq_type = seq_type
        self.seq_side = seq_side
        self.priority = priority
        self.realizations = list()
        self.name = ""
        self.nickname = nickname
        self.update_name()
        return

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Rec_event(" + self.name + ")"

    def add_realization(self, realization):
        """Add a realization to the RecEvent realizations list."""
        self.realizations.append(realization)
        self.update_name()

    def update_name(self):
        """Updates the name of the event (will have no effect if the RecEvent
        has not been modified since the last call).

        """
        if self.event_type == "DinucMarkov":
            self.name = self.event_type + "_" + self.seq_type + "_" + \
                        self.seq_side + "_prio" + \
                        str(self.priority) + "_size" + \
                        str(len(self.realizations) ** 2)
        else:
            self.name = self.event_type + "_" + self.seq_type + "_" + \
                        self.seq_side + "_prio" + \
                        str(self.priority) + "_size" + \
                        str(len(self.realizations))

    def get_realization_vector(self):
        """This methods returns the event realizations sorted by the
        realization index as a list.

        """
        if self.event_type == 'GeneChoice':
            tmp = [""] * len(self.realizations)  # empty(, dtype = str)
        else:
            tmp = numpy.empty(len(self.realizations),
                              dtype=type(self.realizations[0].value))
        # print("Unfinished method get realization vector")
        processed_real_indices = []
        for real in self.realizations:
            if processed_real_indices.count(real.index) == 0:
                if real.name != "":
                    tmp[real.index] = real.name
                else:
                    tmp[real.index] = real.value
                processed_real_indices.append(real.index)
            else:
                print("REALIZATION INDICES ARE DEGENERATE")

        return tmp


class Adjacency_List:
    """Graph adjacency list. Similar to IGoR's C++ Adjacency_List class. Stores
    adjacent events (parents and children) to a given event.

    """
    def __init__(self):
        self.children = list()
        self.parents = list()


class Event_realization:
    """A small class storing for each RecEvent realization its name, value and
    corresponding index.

    """
    def __init__(self, name, value, index):
        self.name = name
        self.value = value
        self.index = index

    def __lt__(self, other):
        return self.index < other.index

    def __gt__(self, other):
        return self.index > other.index


def compute_average_distribution(event_name, model, averaging_list=None):
    """Returns the 1D distribution averaged over the event ancestors.
    Alternatively a list of ancestors can be provided, such that the
    distribution is only averaged over these ancestors.

    """
    for event in model.events:
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
                    if parent not in model.edges:
                        continue
                    parent_parents = model.edges[parent].parents
                    for p in parent_parents:
                        next_parents.append(p)
                        if event_cluster.count(model.get_event(p)) == 0:
                            event_cluster.append(model.get_event(p))
                parents_names = next_parents
                if len(parents_names) == 0:
                    explored = True

            dimension_list = []
            dimension_names_list = []
            for ev in event_cluster:
                dimension_list.append(
                    model.marginals[0][ev.nickname].shape[-1])
                dimension_names_list.append(ev.nickname)

            if event.event_type != "DinucMarkov":

                final_array = numpy.ones(tuple(dimension_list))
                # print(final_array.shape)

                # Reshape all marginals in one big redundant array
                reshaped_marginals = list()

                # Takes care of model 1(non log part)
                for ev in event_cluster:
                    if ev.name not in model.edges:
                        event_parents = []
                    else:
                        event_parents = model.edges[ev.name].parents
                    # print(event_parents)
                    dimension_list = []
                    for e in event_cluster:
                        if (event_parents.count(e.name) > 0) \
                                or (ev.name == e.name):
                            dimension_list.append(
                                model.marginals[0][e.nickname].shape[-1])  # same as getting the number of realizations
                        else:
                            dimension_list.append(1)

                    # print(ev)
                    # print(dimension_list)
                    # print(shape(model1.marginals[ev.nickname]))
                    # print("Dimension list")
                    # print(dimension_list)

                    swapped_marginal = copy.deepcopy(
                        model.marginals[0][ev.nickname])
                    swapped_dimension_names = copy.deepcopy(
                        model.marginals[1][ev.nickname])

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

                for reshaped_marg in reshaped_marginals:
                    final_array *= reshaped_marg

                # if no list given average over all dependencies
                sum_axes = range(0, len(reshaped_marginals[0].shape))
                sum_axes.remove(find(
                    numpy.asarray(dimension_names_list) == event.nickname))
                if averaging_list is not None:
                    # for dependance_name in averaging_list:
                    # sum_axes.remove(find(numpy.asarray(dimension_names_list)==dependance_name))
                    for dimension_name in dimension_names_list:
                        if dimension_name != event.nickname and all(
                                numpy.asarray(
                                        averaging_list) != dimension_name):
                            sum_axes.remove(find(numpy.asarray(
                                dimension_names_list) == dimension_name))
                return final_array.sum(axis=tuple(sum_axes))


# class Model_Marginals:
# Model marginals are for now quite uneasy to study, should be turned into a handy class
# For each event a small object containing the array and another a tuple with the dimension names ordered
def read_marginals_txt(filename, dim_names=False):
    """Reads a model marginals file. Returns a tuple containing a dict
    containing the individual events probabilities indexed by the events
    nicknames and a dict containing the list of dimension names/ordering for
    each event.

    """
    with open(filename, 'r') as f:
        # Model parameters are stored inside a dictionnary of ndarrays
        model_dict = {}
        dimension_names_dict = {}
        element_name = ""
        first = True
        first_dim_line = False
        element_marginal_array = []
        indices_array = []

        for line in f:
            strip_line = line.rstrip('\n')  # Remove end of line character
            if strip_line[0] == '@':
                first_dim_line = True
                if not first:
                    # Add the previous to the dictionnary
                    model_dict[element_name] = element_marginal_array
                else:
                    first = False

                element_name = strip_line[1:]
            # print element_name

            if strip_line[0] == '$':
                # define array dimensions
                coma_index = strip_line.find(',')
                dimensions = []

                # Get rid of $Dim[
                previous_coma_index = 4
                while coma_index != -1:
                    dimensions.append(
                        int(strip_line[previous_coma_index + 1:coma_index]))
                    previous_coma_index = coma_index
                    coma_index = strip_line.find(',', coma_index + 1)

                # Add last dimension and get rid of the closing bracket
                dimensions.append(int(strip_line[previous_coma_index + 1:-1]))

                element_marginal_array = numpy.ndarray(shape=dimensions)

            if strip_line[0] == '#':
                if first_dim_line:
                    dimensions_names = []
                    if len(dimensions) > 1:
                        comma_index = strip_line.find(',')
                        opening_bracket_index = strip_line.find('[')
                        while opening_bracket_index != -1:
                            dimensions_names.append(
                                strip_line[
                                    opening_bracket_index + 1:comma_index])
                            opening_bracket_index = strip_line.find(
                                '[', comma_index)
                            comma_index = strip_line.find(
                                ',', opening_bracket_index)
                    first_dim_line = False
                    dimensions_names.append(element_name)
                    dimension_names_dict[element_name] = dimensions_names

                # update indices
                indices_array = []
                if len(dimensions) > 1:
                    comma_index = strip_line.find(',')
                    closing_brack_index = strip_line.find(']')
                    while closing_brack_index != -1:
                        indices_array.append(int(
                            strip_line[comma_index + 1:closing_brack_index]))
                        opening_bracket_index = strip_line.find(
                            '[', closing_brack_index)
                        comma_index = strip_line.find(
                            ',', opening_bracket_index)
                        closing_brack_index = strip_line.find(
                            ']', closing_brack_index + 1)

            if strip_line[0] == '%':
                # read doubles
                coma_index = strip_line.find(',')
                marginals_values = []

                # Get rid of the %
                previous_coma_index = 0
                while coma_index != -1:
                    marginals_values.append(
                        float(strip_line[previous_coma_index + 1:coma_index]))
                    previous_coma_index = coma_index
                    coma_index = strip_line.find(',', coma_index + 1)

                # Add last dimension and get rid of the closing bracket
                marginals_values.append(
                    float(strip_line[previous_coma_index + 1:]))
                if len(marginals_values) != dimensions[-1]:
                    print("problem")
                element_marginal_array[tuple(indices_array)] = marginals_values
        model_dict[element_name] = element_marginal_array

    return model_dict, dimension_names_dict
