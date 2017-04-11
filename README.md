# **IGoR: Inference and Generation Of Repertoires** #

# *README* #


This repository contains all sources and models useful to infer V(D)J recombination related processes for TCR or BCR sequencing data using **IGoR**

### Quick summary ###

IGoR is a C++ software designed to infer V(D)J recombination related processes from sequencing data such as:

+ Recombination model probability distribution
+ Hypermutation model
+ Best candidates recombination scenarios
+ Generation probabilities of sequences (even hypermutated)

Its heavily object oriented and modular style was designed to ensure long term support and evolvability for new tasks in assessing TCR and BCR receptors features using modern parallel architectures. 

### Version ###
1.0

### Table of Content ###

[TOC]

# Dependencies

+ a C++ compiler supporting OpenMP 3.8 or higher and POSIX Threads (pthread)
+ GSL library : a subpart of the library is shipped with IGoR and will be statically linked to IGoR's executable to avoid dependencies 
+ jemalloc (optional although recommended for full parallel proficiency) memory allocation library: also shipped with IGoR to avoid dependencies issues (requires pthreads)
+ bash

# Install
IGoR uses the autotools suite for compilation and installation in order to ensure portability to many systems. 

## Linux
Widely tested on several Debian related distros.
Install gcc/g++ if not already installed (although can also be compiled using icc for instance).
With the command line go to IGoR's root directory and simply type `./configure`. This will make various check on your system and create makefiles compatible with your system configuration. Once over, type `make` to compile the sources and obtain IGoR's executable.

## MacOS
MacOS is shipped with another compiler (Clang) when installing Xcode that is called upon calling gcc and is not supporting OpenMP. In order to use gcc and compile with it an OpenMP application you will first need to download Macports and install gcc from there.

Once done, as for Linux simply go to IGoR's root folder and type `./configure;make` 

## Windows (not tested)
The configure script relies on bash to work. A first step is to download a bash interpreter (such as Cygwin or MinGW). Open the command line of the one of your choice and use `./configure;make`

# Workflow
As a preprocessing step IGoR first needs to align the genomic templates to the read (`-align`) before exploring all putative recombination scenarios for this read. 
After aligning IGoR can be used to infer a recombination model (`-infer`), evaluate sequences statistics (`-evaluate`) using an already inferred model.
Synthetic sequences can be generated from a learned model (as one supplied by IGoR, or one inferred de novo through the `-infer` command) with the `-generate` command. 

# Command line tools
Although the full flexibility of IGoR is reachable through C++ highlevel functions (*see next section*) we provide some command line options to perform most frequent tasks on immune receptor sequences.

Command options are nested arguments, the general organization of the commands follows `-arg1 --subarg1 ---subsubarg1` to reach the different levels.

## General

### Using predefined genomic templates and models
IGoR is shipped with a set of genomic templates and already inferred models  from *cite IGoR paper, and hopefully more to come *.

** Some installations issues still need to be figured out, in order to use the predefined models please execute IGoR in the folder in which the executable was compiled.**
 
Available options are listed below:

| Species | Chains |
| :------------- | :------------------------------ | 
| human | alpha,beta,heavy |

If you are working on datasets not present in this list and would kindly agree to contribute to this database please contact us.

### General commands summary

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `-set_wd /path/to/dir/`      | Sets the working directory to */path/to/dir/*, default is */tmp*. ** This should be an already existing directory and will not be created by IGoR **   |
| `-threads N`   | Sets the number of OpenMP threads to *N* for alignments and inference     |
| `-stdout_f /path/to/file`  | Redirects the standard output to the file */path/to/file*  |
| `-read_seqs /path/to/file`  | Reads the input sequences file */path/to/file* and reformat it in the working directory. **This step is necessary for running any action on sequences using the command line**. Can be a fasta file or a text file with one sequence per line (format recognition is based on the file extension). Providing this file will create a semicolon separated file with indexed sequences in the *align* folder.|
|`-batch batchname`| Sets the batch name. This name will be used as a prefix to alignment/indexed sequences files, output, infer, evaluate and generate folders.|
| `-chain chainname` | Selects a model and a set of genomic template according to the value. Possible values for `chainname` are: `alpha`, `beta`, `light`, `heavy_naive`, and `heavy_memory`. **This needs to be set in order to use provided genomic templates/model**|
| `-species speciesname`| Selects a species from the set of predefined species. Possible values are: `human`.**This needs to be set in order to use provided genomic templates/model** |
|`-set_genomic --*gene* /path/to/file.fasta`| Set a set of custom genomic templates for gene *gene* (possible values are V,D and J) with a list of genomic templates contained in the file */path/to/file.fasta* in fasta format. ** When using this option you will need to re-infer a model since the genomic templates will no longer correspond to the ones contained in the reference models **|
| `-set_custom_model /path/to/model_parms.txt /path/to/model_marginals.txt` | Use a custom model as a baseline for inference or evaluation. **Note that this will override  custom genomic templates for inference and evaluation**|
|`-load_last_inferred`| Using this command will load the last inferred model (folder *inference/final_xx.txt*) as a basis for a new inference, evaluation or generation of synthetic sequences |
| `-run_demo`  |  Runs the demo code on 300 sequences of 60bp TCRs (mostly a sanity run check) |
| `-run_custom`  |  Runs the code inside the custom section of the main.cpp file |

### Working directory
This is where all IGoR outputs will appear. Specific folders will be created for alignments, inference, evaluation and outputs.

## Alignments

### Algorithm

Performs Smith-Waterman alignments of the genomic templates. Using a slight alteration of the Smith-Waterman score matrix, we enforce that V can only be deleted on the 3' side and J on the 5' side (thus enforcing the alignment on the other side until the end of the read or of the genomic template). D is aligned using a classical Smith-Waterman local alignment approach allowing gene deletions on both sides.

### Alignment commands summary
 
Alignment of the sequences is performed upon detection of the `-align` switch in the command line. For each gene, alignment parameters can be set using `--V`,`--D` or `--J`. **Specifying any of those three argument will cause to align only the specified genes**. In order to specify a set of parameters for all genes or force to align all genes the argument `--all` should be passed. The arguments for setting the different parameters are given in the table below.

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `---thresh X`  | Sets the score threshold for the considered gene alignments to *X*. Default is 50.0 for V, 15.0 for D and 15.0 for J |
| `---matrix path/to/file` | Sets the substitution matrix to the one given in the file. Must be ZZZ delimited. Default is a NUC44 matrix with stronger penalty on errors (5,-14) ** Not fully supported yet **|
| `---gap_penalty X` | Sets the alignment gap penalty to X. Default is 50.0 |
| `---best_only` | If *true* only keep the best alignment per gene/allele. If *false* outputs all alignments above the score threshold. Default is *true* for V and J, and *false* for D. |
| `---offset_bounds M N` | Constrains the possible positions of the alignments. The offset is defined as the position on the read to which the first nucleotide of the genomic template aligns (can be negative, e.g for V for which most of the V is on the 5' of the read and cannot be seen)  |

### Alignment output files summary
Upon alignment the alignment parameters/dates/filenames will appended to the *aligns/aligns_info.out* file for easy traceability.

Alignment files are semicolon separated files. For each alignment of a genomic template to a sequence the following fields are given:

| Field | Description |
| :------------- | :------------------------------ |
| seq_index | The sequence index the alignment corresponds to in the *indexed_sequences.csv* file. |
| gene_name | The gene name as provided in the genomic template file |
| score | SW alignment score |
| offset | The index of the first letter of the (undeleted) genomic template on the read as described in the previous section. |
| insertions | Indices of the alignment inserted nucleotides (relative to the read) |
| deletions | Indices of the alignment deleted nucleotides (relative to the genomic template) |
| mismatches | Indices of the alignment mismatches (relative to the read) |
| length | Length of the SW alignment (including insertions and deletions) |
| 5_p_align_offset | Offset of the first nucleotide of the SW alignment (relative to the read) |
| 3_p_align_offset | Offset of the last nucleotide of the SW alignemnt (relative to the read) |


## Inference and evaluation

### Inference and evaluation commands
The inference is reached using the command `-infer`. Logs and models parameters values for each iteration will be created in the folder *inference* of the working directory (or *batchname_inference* if a batchname was supplied). 

Sequence evaluation is reached using the command `-evaluate`. This is the same as performing an iteration of the Expectation-Maximization on the whole dataset and thus accepts the same arguments as `-infer` for arguments related to the precision of the algorithm. The logs of the sequences evaluation are created in the folder *evaluate* (or *batchname_evaluate* if a batchname was supplied).

** Note that -infer and -evaluate are mutually exclusive in the same command since it brings ambiguity reagarding which model should be used for each **

Optional parameters are the following:

| Command line argument | Description                    | Available for |
| :------------- | :------------------------------ | :-----------|
| `--N_iter N`  | Sets the number of EM iterations for the inference to N| inference |
| `--L_thresh X`  | Sets the sequence likelihood threshold to X. | inference & evaluation |
| `--P_ratio_thresh X`  | Sets the probability ratio threshold to X. This influences how much the tree of scenarios is pruned. Setting it 0.0 means exploring every possible scenario (exact but very slow), while setting it to 1.0 only explores scenarios that are more likely than the best scenario explored so far (very fast but inaccurate). This sets a trade off between speed and accuracy, the best value is the largest one for which the likelihood of the sequences almost doesn't change when decreasing it further.  | inference & evaluation |
| `--MLSO`  | Runs the algorithm in a 'Viterbi like' fashion. Accounts for the Most Likely Scenario Only (as fast as using a probability ratio threshold of 1.0) | inference & evaluation |
|`--infer_only eventnickname1 eventnickname2`| During the inference only the the parameters of the events with nicknames listed are updated | inference |
|`--not_infer eventnickname1 eventnickname2`| Opposite command to the one above, will fix the parameters of the listed events | inference |
|`--fix_err`| In the same vein as the two commands above, this one will fix the parameters related to the error rate. | inference |
| `--subsample N` | Perform the action on a random subsample of *N* sequences. ** Obviously N should be < to the total number of sequences available ** | inference & evaluation |

### Inference and evaluation output
Upon inferring or evaluating several files will be created in the corresponding folder.

1.*\*_parms.txt* files contain information to create Model_Parms C++ objects. It encapsulates information on the model events, realizations and topology. 
2.*\*_marginals.txt* files contain information to create Model_Marginals C++ objects. It encapsulate parameters of the probability distribution underlying the recombination model.
3.*inference_info.out* contains the inference parameters for traceability
4.*inference_logs.txt* contains some information on each sequence for each iteration. This is a useful tool to debug inference troubleshoots. 
5.*likelihoods.out* contains the likelihood information for a given dataset.

### Inference and evaluation Troubleshoots
Although the inference/evaluation generally run smoothly we try to list out some possible troubleshoots and corresponding solutions.

| Issue | Putative solution                    | 
| :------------- | :------------------------------ | 
|map_base::at() exception | This exception is most likely thrown by a Gene_Choice event in the inference. Try/Catch handling is runtime costly thus some checks are not performed on the fly. Explanation: This is most likely the inference receiving a genomic template whose name does not exist in the model realizations. Solution: make sure the genomic templates (and their names) used for alignments correspond to those contained in your model file. |
| All 0 output | All marginal files contains 0 parameters after one iteration. All sequences have zero likelihood in the *inference_logs.txt* file. Explanation: none of the scenarios had a sufficiently high likelihood to reach the likelihood threshold. Solution: use the `--L_thresh` argument to decrease the likelihood threshold, if the code becomes utterly slow see below. ** In general while inferring one should make sure not too many sequences are assigned a zero likelihood since it would introduce a systematic bias in the learned distribution ** |
| Extreme slowness | Runtimes are very far from the ones given in *cite IGoR*. Check the mean number of errors in the *inference_logs.txt* file. If these numbers are higher than you would expect from your data (e.g if you are not studying hypermutated data) check your alignments statistics. A possible explanation would be an incorrect setting of the alignment offsets bounds |


## Outputs 
Outputs are scenario/sequence statistics each presented below. They are all written in the *output* folder.
The different outputs are detailed in the next sections.

### Best scenarios
*Output the N best scenarios for each sequence*

Use command `--scenarios N`

### Generation probability
*Estimates the probability of generation of the error free/unmutated ancestor sequence*

Use command `--Pgen`

### Coverage
*Counts for each genomic nucleotide how many times it has been seen and how many times it was mutated/erroneous*

Use command `--coverage`
## Sequence generation
Reached using the command `-generate N` where *N* is the number of sequences to be generated. The number of sequences to generate must be passed before optional arguments. Optional parameters are the following:

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `--noerr`  | Generate sequences without sequencing error (the rate and the way those errors are generated is controlled by the model error rate)|
| `--name myname`  | Prefix for the generated sequences filenames. ** Note that setting the *batchname* will change the generated sequences folder name, while setting *--name* will change the file names. ** |
| `--seed X`  | Impose *X* as a seed for the random sequence generator. By default a random seed is obtained from the system. |

## Command examples
Here we give a few command examples for a typical workflow.


# C++
Although a few command line options are supplied for basic use of IGoR, its full modularity can be used through high level C++ functions on which all previous command lines are built. A section of the main.cpp file is dedicated to accept user supplied code and can be executed using `-run_custom` command line when launching IGoR from the shell. An example of the workflow is given in the *run demo* section and the full Doxygen generated documentation is available as PDF. For any question please contact us.

Good practice would be to append the C++ code in the main in the scope where "//Write your custom procedure here" is written. This part of the code is reachable using the `-run_custom` command line argument. This is done so that even after appending some custom code the command line interface is still usable.

# Python
A set of Python modules are shipped with Igor in order to parse IGoR's outputs (alignments,models etc)
For further versions a Python/Cython interface for IGoR might be supplied 

# Contribute

* Your feedbacks are valuable, please send your comments about usability and new features you would like to see  
* Code contribution: IGoR was designed to be modular and evolve, please get in touch if you would like to do something new with your data and would like some more guidance on the code structure


# Contact 

For any question please file an issue on github or email quentin.marcou@lpt.ens.fr