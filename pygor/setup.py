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

from setuptools import setup

with open('LICENSE') as f:
    license = f.read()

setup(name='pygor',
      version='1.2.0',
      description='Module for parsing IGoR output (alignments, models etc).',
      url='https://github.com/qmarcou/IGoR',
      author='Quentin Marcou',
      author_email='quentin.marcou@lpt.ens.fr',
      license=license,
      packages=[
          'pygor',
          'pygor.aligns',
          'pygor.counters',
          'pygor.counters.bestscenarios',
          'pygor.counters.coverage',
          'pygor.models',
          'pygor.utils'
      ],
      install_requires=[
          'pandas>=0.20.3',
          'numpy>=1.13.1',
          'scipy>=1.0.0',
          'matplotlib>=2.0.2',
          'biopython>=1.70'
      ],
      zip_safe=False)
