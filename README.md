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

# Dependencies

+ GCC (sufficiently recent to work with C++11 standards and OpenMP 3.8)
+ GSL library : currently working on shipping it with IGoR
+ jemalloc (optional although recommended for full parallel proficiency)

# Install

## Linux
Widely tested on several Debian related distros.
Install gcc/g++ if not already installed.

## MacOS
MacOS is shipped with another compiler (Clang) when installing Xcode that is called upon calling gcc and is not supporting OpenMP. In order to use gcc and compile with it an OpenMP application you will first need to download Macports and install gcc from there.


* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

# Workflow
As a preprocessing step IGoR first needs some alignments of the genomic templates to the read before exploring all putative recombination scenarios for this read.

# Command line tools
Although the full flexibility of IGoR is reachable through C++ highlevel functions (*see next section*) we provide some command line options to perform most frequent tasks on immune receptor sequences.

## General
| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `-set_wd /path/to/dir/`      | Sets the working directory to */path/to/dir/*, default is /tmp   |
| `-threads N`   | Sets the number of OpenMP threads to *N* for alignments and inference     |
| `-stdout_f /path/to/file`  | Redirects the standard output to the file */path/to/file*  |
| `-read_seqs /path/to/file`  | Reads the input sequences file */path/to/file* and reformat it in the working directory. **This step is necessary for running any action on sequences using the command line**. Can be a fasta file or a text file with one sequence per line (format recognition is based on the file extension). |
| `-chain --chainname` | Selects a model and a set of genomic template according to the value. Possible values for `--chainname` are: `--alpha`, `--beta`, `--light`, `--heavy_naive`, and `--heavy_memory`. **This needs to be set in order to use provided genomic templates/model**
| `-run_demo`  |  Runs the demo code on 300 sequences of 60bp TCRs (mostly a sanity run check) |

### Working directory
This is where all IGoR outputs will appear. Specific folders will be created for alignments, inference , evaluation and outputs.

## Alignments
Performs Smith-Waterman alignments of the genomic templates. Using a slight alteration of the smith waterman score matrix, we enforce that V can only be deleted on the 3' side and J on the 5' side (thus enforcing the alignment on the other side until the end of the read or of the genomic template). D is aligned using a classical smith waterman local alignment approach allowing gene deletions on both sides.
Alignment of the sequences is performed upon detection of the `-align` switch in the command line. For each gene, alignment parameters can be set using `--V`,`--D` or `--J`. The arguments for setting the different parameters are given in the table below.

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `---thresh X`  | Sets the score threshold for the considered gene alignments to *X*. Default is ZZ for V, ZZ for D and ZZ for J |
| `---matrix path/to/file` | Sets the substitution matrix to the one given in the file. Must be ZZZ delimited. Default is a NUC44 matrix with stronger penalty on errors (5,-14) |
| `---gap_penalty X` | Sets the gap penalty to X |
| `---best_only` | |
| `---offset_bounds M N` | Constrains the possible positions of the alignments. The offset is defined as the position on the read to which the first nucleotide of the genomic template aligns (can be negative, e.g for V for which most of the V is on the 5' of the read and cannot be seen)  |

## Inference
The inference is reached using the command `-infer`. Logs and models parameters values for each iteration will be created in the folder *inference* of the working directory. Optional parameters are the following:

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `--N_iter N`  | Sets the number of EM iterations for the inference to N|
| `--L_thresh X`  | Sets the sequence likelihood threshold to X. |
| `--P_ratio_thresh X`  | Sets the probability ratio threshold to X. This influences how much the tree of scenarios is pruned. Setting it 0.0 means exploring every possible scenario (exact but very slow), while setting it to 1.0 only explores scenarios that are more likely than the best scenario explored so far (very fast but inaccurate). This sets a trade off between speed and accuracy, the best value is the largest one for which the likelihood of the sequences almost doesn't change when decreasing it further.  |
| `--viterbi`  | Runs the algorithm in a 'Viterbi like' fashion. Only accounts for the most likely scenario (as fast as using a probability ratio threshold of 1.0) |

### Troubleshoots
map base at exception => check genomic templates (explain try catch expensive)
run smoothly but all 0: alignments!!

## Evaluate
Reached using the command `-evaluate`. This is the same as performing an iteration of the inference on the whole dataset and thus accepts the same arguments as `-infer` except for `--N_iter`. The logs of the sequences evaluation are created in the folder *eval*.

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
Reached using the command `-generate N` where *N* is the number of sequences to be generated. Optional parameters are the following:

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `--noerr`  | Generate sequences without sequencing error (the rate and the way those errors are generated is controlled by the model error rate)|
| `--name batchname`  | Batch name for the generated sequences filename (useful if one needs to generate several sets of generated sequences). |


# C++
Although a few command line options are supplied for basic use of IGoR, its full modularity can be used through high level C++ functions on which all previous command lines are built. A section of the main.cpp file is dedicated to accept user supplied code and can be executed using `-custom` command line when launching IGoR from the shell. An example of the workflow is given in the *run demo* section and the full Doxygen generated documentation is available as PDF. For any question please contact us.

# Python
A set of Python modules are shipped with Igor in order to parse IGoR's outputs (alignments,models etc)
For further versions a Python/Cython interface for IGoR might be supplied 

# Contribute

* Your feedbacks are valuable, please send your comments about usability and new features you would like to see  
* Code contribution: IGoR was designed to be modular and evolve, please get in touch if you would like to do something new with your data and would like some more guidance on the code structure


# Contact 

For any question please file an issue on github or email quentin.marcou@lpt.ens.fr