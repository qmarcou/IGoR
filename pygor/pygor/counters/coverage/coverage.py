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

import numpy
import pandas

from ...utils.utils import get_str_asarray


def compute_mutation_frequency(x, threshold):
    """Computes mutation frequency for positions with coverage larger than
    threshold.

    """
    fraction_list = []
    for cov, err in zip(x.coverage, x.errors):
        if cov > threshold:
            fraction_list.append(err / cov)
        else:
            fraction_list.append(numpy.NaN)
    return fraction_list


# Read gene coverage and errors data
def read_coverage_and_errors_file(filename, get_diag_of_N_dim=None):
    """Reads coverage and error counter file as pandas.DataFrame. For
    Ndimensional joint coverage and errors the diagonal(equivalent to 1D) can
    be extracted by passing the number of dimensions as an argument.

    """
    raw_read = pandas.read_csv(filename, sep=';')

    # Convert coverage strings into arrays
    tmp_cov = raw_read.apply(lambda x: get_str_asarray(x.coverage), axis=1)
    # Convert errors strings into arrays
    tmp_err = raw_read.apply(lambda x: get_str_asarray(x.errors), axis=1)

    if get_diag_of_N_dim is not None:
        for i in range(0, len(tmp_cov)):
            single_dim_len = int(float(len(tmp_cov.iloc[i])) ** (
                        1.0 / float(get_diag_of_N_dim)))
            reshape_tuple = tuple([single_dim_len] * get_diag_of_N_dim)
            # print(reshape_tuple)
            # print(len(tmp_cov))

            # Get the diagonal
            tmp_cov.iloc[i] = list(numpy.asarray(tmp_cov.iloc[i]).reshape(
                reshape_tuple).diagonal())
            tmp_err.iloc[i] = list(numpy.asarray(tmp_err.iloc[i]).reshape(
                reshape_tuple).diagonal())

    raw_read.coverage = tmp_cov
    raw_read.errors = tmp_err

    return raw_read
