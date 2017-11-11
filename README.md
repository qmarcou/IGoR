# **IGoR: Inference and Generation Of Repertoires** #

# *README* #


This repository contains all sources and models useful to infer V(D)J recombination related processes for TCR or BCR sequencing data using **IGoR**

### Quick summary ###

IGoR is a C++ software designed to infer V(D)J recombination related processes from sequencing data such as:

+ Recombination model probability distribution
+ Hypermutation model
+ Best candidates recombination scenarios
+ Generation probabilities of sequences (even hypermutated)

The following paper describes the methodology, performance tests and some new biological results obtained with IGoR:

[*IGoR: A Tool For High-Throughput Immune Repertoire Analysis.*][igor_bioarxiv] (2017) Quentin Marcou, Thierry Mora, Aleksandra M. Walczak

[igor_bioarxiv]: https://arxiv.org/abs/1705.08246

Its heavily object oriented and modular style was designed to ensure long term support and evolvability for new tasks in assessing TCR and BCR receptors features using modern parallel architectures. 

### Version ###
Latest released version: 1.1.0

### Table of Content ###

[TOC]

# Dependencies

+ a C++ compiler supporting OpenMP 3.8 or higher and POSIX Threads (pthread) such as GCC
+ [GSL library](https://www.gnu.org/software/gsl/) : a subpart of the library is shipped with IGoR and will be statically linked to IGoR's executable to avoid dependencies 
+ [jemalloc](http://jemalloc.net/) (optional although recommended for full parallel proficiency) memory allocation library: also shipped with IGoR to avoid dependencies issues (requires a pthreads compatible compiler)
+ bash
+ autotools suite if building from unpackaged sources

# Install
IGoR uses the autotools suite for compilation and installation in order to ensure portability to many systems. 

*Installing from package releases (recommended)*

First download the latest released package on the download page (on the left). Extract the files from the archive. 

*Installing from unpackaged sources*

For this you will have to get git, the autotools suite, pandoc, doxygen and the latex suite software installed. Note that this is the most convenient way to keep IGoR up-to-date but involves a bit more installation steps.
Using *git*, clone the repository where you desire. Go in the created directory and run the *autogen.sh* bash script. This will create the *configure* script. Upon this stage the installation rules are the same as for packaged developper sources.
From *git* you can chose among two branches: the *master* branch corresponds to the latest stable (packaged) release, the dev branch is the most up to date branch including current developpments until they are issued in the next release. The *dev* branch is therefore more bug prone, however this is the natural branch for people ready to help with developpment (even only by functionality testing).

A (sadly) non exhaustive list of potential installation troubleshoots follows in the next section. If your problem is not referenced there please [contact](<quentin.marcou@lpt.ens.fr> "myadress") us. If you end up finding a solution by yourself please help us append it to the following list and help the user community.

## Linux
Widely tested on several Debian related distros.
Install gcc/g++ if not already installed (note that another compiler could be used).
With the command line go to IGoR's root directory and simply type `./configure`. This will make various check on your system and create makefiles compatible with your system configuration. Many options can be appended to ./configure such as `./configure CC=gcc CXX=g++` to enforce the use of gcc as compiler. Once over, type `make` to compile the sources (this will take a few minutes). **IGoR's executable will appear in the igor_src folder**

Finally in order to access all IGoR's features, install IGoR by typing `make install`. This will install IGoR's executable, supplied models and manual on your system. If you do not have administrator privileges, IGoR can be installed locally in the folder of your choice by passing **--prefix=/your/custom/path** upon calling the configure script (e.g `./configure --prefix=$HOME`). Other configure options can be accessed using ./configure -h. 

## MacOS
MacOS is shipped with another compiler (Clang) when installing Xcode that is called upon calling gcc (through name aliasing) and is not supporting OpenMP. In order to use gcc and compile with it an OpenMP application you will first need to download Macports or Homebrew and install gcc from there.

First if not already present on your system install XCode through the application store. 

Macports can be found [here][macports_site]. Download and install the version corresponding to your MacOS version.
[macports_site]: https://www.macports.org/install.php

Once installed, use Macports to install GCC:
```
sudo port selfupdate #Update macports database
sudo port install gcc6 #install gcc version 6
```
The full list of available GCC versions is available [here][macports_gccs], select a sufficiently recent one to get C++11 standards enabled. In order to set GCC as your default compiler use the following commands:
[macports_gccs]: https://www.macports.org/ports.php?by=name&substr=gcc

```
port select --list gcc #Will list the versions of gcc available on your system
sudo port select --set gcc mp-gcc6 #set the one you wish to have as default call upon using the gcc command
```

If you prefer to use Homebrew over Macports, it can be downloaded and installed [here](https://brew.sh/).

Then install GCC using the following command:

`brew install gcc`

**Note: if you decide to use Homebrew you should apparently refrain yourself from assigning the newly installed gcc to the `gcc` command(see [this page](http://docs.brew.sh/Custom-GCC-and-cross-compilers.html) for more details). You will thus have to pass the correct compiler instructions to the configure script with the *CC* and *CXX* flags.** 

Once done, as for Linux simply go to IGoR's root folder and type `./configure` to create the appropriate makefiles and compile using `make` (this will take a few minutes). Finally install IGoR using the `make install` command.


## Windows (not tested)
The configure script relies on bash to work. A first step would be to download a bash interpreter (such as Cygwin or MinGW) and a compiler. Open the command line of the one of your choice and use `./configure;make`

## Troubleshoots 
Here is a list of some install troubleshoots that have been reported and their corresponding solution

| Issue | Reason | Solution |
| :------------- | :------------------------------ | :------------------------------ |
|In file included from Aligner.cpp:8: /n ./Aligner.h:19:10: fatal error: 'omp.h' file not found /n #include <omp.h>| The compiler used is not supporting OpenMP | Make sure you have an OpenMP compatible compiler installed (such as GCC). If such a compiler is installed make sure the right compiler is called upon compiling. In order to specify a specific compiler to use (such as mc-gcc6 for macport installed gcc under MacOS) pass the following option upon executing the configure script: `./configure CC=mc-gcc6 CXX=mc-g++6`. The *CC* option will specify the C compiler to use to compile jemalloc and gsl, while *CXX* specifies the C++ compiler to use to compile IGoR sources. |
|  *aclocal-1.15: command not found*; *WARNING: 'aclocal-1.15' is missing on your system.*; *make: *** [aclocal.m4] Error 127* | The *configure* script relies on file timestamps to assess whether it is up to date. These time stamps might be compromised when extracting files from the archive. | Run the following command in IGoR root directory: `touch configure.ac aclocal.m4 configure Makefile.* */Makefile.* */*/Makefile.*` |
| *.libs/sasum.o: No such file or directory* error at compile time | Unknown | Running `make clean;make` will fix this issue |
|*undefined reference to symbol 'clock_gettime@@GLIBC_2.2.5'* at link time| Jemalloc used an extra library to extract system time | Run the last command printed to the screen (*g++ -std=gnu++11 -I./../libs/jemalloc/include/jemalloc -I./../libs/gsl_sub -fopenmp ...... -lpthread -ldl -fopenmp*) and add -lrt after -ldl. This will be automated and fixed soon |
| *src/jemalloc.c:241:1: error: initializer element is not constant* ; *static malloc_mutex_t init_lock = MALLOC_MUTEX_INITIALIZER;* | Might be related to MacOS Sierra? | Unknown |

# Workflow
As a preprocessing step IGoR first needs to align the genomic templates to the read (`-align`) before exploring all putative recombination scenarios for this read. 
After aligning IGoR can be used to infer a recombination model (`-infer`), evaluate sequences statistics (`-evaluate`) using an already inferred model.
Synthetic sequences can be generated from a learned model (as one supplied by IGoR, or one inferred de novo through the `-infer` command) with the `-generate` command. 

# Command line tools
Although the full flexibility of IGoR is reachable through C++ highlevel functions (*see next section*) we provide some command line options to perform most frequent tasks on immune receptor sequences.

Command options are nested arguments, the general organization of the commands follows `-arg1 --subarg1 ---subsubarg1` to reach the different levels.

## General

### Using predefined genomic templates and models
IGoR is shipped with a set of genomic templates and already inferred models from [[1][igor_bioarxiv]].

**In order to use the predefined models and demo IGoR must have been installed on your system.**
 
Available options are listed below:

| Species | Chains |
| :------------- | :------------------------------ | 
| human | alpha,beta,heavy |

If you are working on datasets not present in this list refer to the *Advanced usage section* and/or contact us for assistance. Help us filling this database for other users and share the resulting models with us!

### General commands summary

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `-h` or `-help`   | Displays IGoR's manual. Alternatively one could use `man igor`.  |
| `-v` or `-version`      |  Displays IGoR's installed version number. |
| `-set_wd /path/to/dir/`      | Sets the working directory to */path/to/dir/*, default is */tmp*. ** This should be an already existing directory and will not be created by IGoR **   |
| `-threads N`   | Sets the number of OpenMP threads to *N* for alignments and inference/evaluation. By default IGoR will use the maximum number of threads.     |
| `-stdout_f /path/to/file`  | Redirects the standard output to the file */path/to/file*  |
| `-read_seqs /path/to/file`  | Reads the input sequences file */path/to/file* and reformat it in the working directory. **This step is necessary for running any action on sequences using the command line**. Can be a fasta file, a csv file (with the sequence index as first column and the sequence in the second separated by a semicolon ';') or a text file with one sequence per line (format recognition is based on the file extension). Providing this file will create a semicolon separated file with indexed sequences in the *align* folder.|
| `-batch batchname`| Sets the batch name. This name will be used as a prefix to alignment/indexed sequences files, output, infer, evaluate and generate folders.|
| `-chain chainname` | Selects a model and a set of genomic template according to the value. Possible values for `chainname` are: `alpha`, `beta`, `light`, `heavy_naive`, and `heavy_memory`. **This needs to be set in order to use provided genomic templates/model**|
| `-species speciesname`| Selects a species from the set of predefined species. Possible values are: `human`.**This needs to be set in order to use provided genomic templates/model** |
| `-set_genomic --*gene* /path/to/file.fasta`| Set a set of custom genomic templates for gene *gene* (possible values are --V,--D and --J) with a list of genomic templates contained in the file */path/to/file.fasta* in fasta format. If the set of provided genomic templates is already fully contained (same name and same sequence) in the loaded model (default, custom, last_inferred), the missing ones will be set to zero probability keeping the ratios of the others. For instance providing only one already known genomic template will result in a model with the considered gene usage to be 1.0, all others set to 0.0. **When using this option and introducing new/modified genomic templates, the user will need to re-infer a model since the genomic templates will no longer correspond to the ones contained in the reference models, the model parameters are thus automatically reset to a uniform distribution.** |
| `-set_CDR3_anchors --*gene*` | Load a CSV file containing the index of the CDR3 anchors for the *gene*(--V or --J). The index should correspond to the first letter of the cystein(for V) or tryptophane/phenylalanin (for J) for the nucleotide sequence of the gene. |
| `-set_custom_model /path/to/model_parms.txt /path/to/model_marginals.txt` | Use a custom model as a baseline for inference or evaluation. **Note that this will override  custom genomic templates for inference and evaluation**. Alternatively, providing only the model parameters file will lead IGoR to create model maginals initialized to a uniform distribution. |
| `-load_last_inferred`| Using this command will load the last inferred model (folder *inference/final_xx.txt*) as a basis for a new inference, evaluation or generation of synthetic sequences |
| `-run_demo`  |  Runs the demo code on 300 sequences of 60bp TCRs (mostly a sanity run check) |
| `-run_custom` |  Runs the code inside the custom section of the main.cpp file |
| `-subsample N` | Perform actions on a random subsample of *N* sequences. **This flag will have different effects depending on the supplied commands:** if the `-read_seqs` command is used, the resulting indexed sequence file will be a subsample of sequences contained in the original file. Else, if the `-align` command is used the alignments will be performed on a subsample of the indexed sequences. Else, if the `-evaluate` or `-infer` command is used the inference will be run on a subsample of the indexed sequences. *Obviously N should be < to the total number of sequences available. The `-subsample` flag should be used in only one command of a pipeline, see the Command example section for details.* | 

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
| `---matrix path/to/file` | Sets the substitution matrix to the one given in the file. Must be *\',\'* delimited. Default is a NUC44 matrix with stronger penalty on errors (5,-14) |
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
|`--infer_only eventnickname1 eventnickname2`| During the inference only the parameters of the events with nicknames listed will be updated | inference |
|`--not_infer eventnickname1 eventnickname2`| Opposite command to the one above, will fix the parameters of the listed events | inference |
|`--fix_err`| In the same vein as the two commands above, this one will fix the parameters related to the error rate. | inference |


### Inference and evaluation output
Upon inferring or evaluating several files will be created in the corresponding folder.

#### Model parameters files
*\*_parms.txt* files contain information to create Model_Parms C++ objects.
It encapsulates information on the individual model events, their possible realizations, the model's graph structure encoding events conditional dependences and the error model information. All fields are semi colon separated.
The different sections of the files are delimited by an *\"@\"* symbol, each further subdivided as follows:

- *\"@Event_list\"* introduces the section in which the recombination events (i.e the Bayesian Network/graph nodes) are defined. 
	- *\"#\"* introduces a new recombination event (or node). The line contains 4 fields: 
		- the event type (*GeneChoice*, *Deletion*, *Insertion*, *DinucMarkov*)
		- the targeted genes (*V_gene*, *VD_genes*, *D_gene*, *DJ_genes*, *J_gene*, *VJ_genes*)
		- the gene side (*Five_prime*, *Three_prime*, *Undefined_side*)
		- the event priority: an integer influencing the order in which events are processed during the inference such that events with high priority are preferentially processed earlier.
		- the event nickname
	- *\"%\"* introduces a new event realization. Depending on the recombination event, the first fields will define the realization name and/or values (e.g gene name and gene sequence for *GeneChoice* or number of deletions for *Deletion*) while the final field always denotes the realization's index on the probability array. **This index is automatically assigned by IGoR upon addition of an event realization, changing it will cause undefined behavior.** See the *Advanced usage* section of this README for more information on how to add/remove events realizations.
- \"@Edges\" introduces the section in which the conditionnal dependencies (i.e graph directed edges) are defined.
	- *\"%parent;child\"* introduces a new directed edge/conditional dependence between the parent and child event.
- \"@ErrorRate\" introduces the section in which the error model is defined.
	- *\"#\"* introduces a new error model, the first field defining the error model type and subsequent fields other meta parameters of the error model 
		- *\"%\"* introduces the parameters values linked to the actual error/mutation rate.



#### Model marginals files
*\*_marginals.txt* files contain information to create Model_Marginals C++ objects.
It encapsulates the probabilities for each recombination event's realization. 
As for the model parameters files, the marginals files are are sectionned by special characters as follows:

- *\"@\"* introduces the recombination event's nickname the following probabilities are referring to.
- *\"$Dim\"* introduces the dimensions of the event and its conditional dimensions probability array. By convention the last dimension refers to the considered event dimension.
- *\"#\" introduces the indices of the realizations of the parent events and their nickname corresponding to the following 1D probability array
- *\"%\" introduces the 1D probability array for all of the considered event realizations for fixed realizations of the parents events whose indices were given in the previous line.

Python functions are provided to read such files along with the corresponding model parameters file within the GenModel object.

#### Inference information file
*inference_info.out* contains the inference parameters/date/time for traceability and potential error messages.

#### Inference logs file
*inference_logs.txt* contains some information on each sequence for each iteration. This is a useful tool to debug inference troubleshoots. 

#### Model likelihood file
*likelihoods.out* contains the likelihood information for a given dataset.

### Inference and evaluation Troubleshoots
Although the inference/evaluation generally run smoothly we try to list out some possible troubleshoots and corresponding solutions.

| Issue | Putative solution                    | 
| :------------- | :------------------------------ | 
|map_base::at() exception | This exception is most likely thrown by a Gene_Choice event in the inference. Try/Catch handling is runtime costly thus some checks are not performed on the fly. Explanation: This is most likely the inference receiving a genomic template whose name does not exist in the model realizations. Solution: make sure the genomic templates (and their names) used for alignments correspond to those contained in your model file. |
| All 0 output | All marginal files contains 0 parameters after one iteration. All sequences have zero likelihood in the *inference_logs.txt* file. Explanation: none of the scenarios had a sufficiently high likelihood to reach the likelihood threshold. Solution: use the `--L_thresh` argument to decrease the likelihood threshold, if the code becomes utterly slow see below. ** In general while inferring one should make sure not too many sequences are assigned a zero likelihood since it would introduce a systematic bias in the learned distribution ** |
| Extreme slowness | Runtimes are very far from the ones given in [the original article][igor_bioarxiv]. Check the mean number of errors in the *inference_logs.txt* file. If these numbers are higher than you would expect from your data (e.g if you are not studying hypermutated data) check your alignments statistics. A possible explanation would be an incorrect setting of the alignment offsets bounds |


## Outputs 
Outputs or Counters in the C++ interface are scenario/sequence statistics, each individually presented below. They are all written in the *output* folder (or *batchname_output* if a batchname was supplied).

In order to specify outputs use the `-output` argument, and detail the desired list of outputs. Outputs are tied to the exploration of scenarios and thus require to have `-infer` or `-evaluate` in the same command. Note that although it might be interesting to track some outputs during the inference for debugging purposes, best practice would be to use it along with evaluation. 

The different outputs are detailed in the next sections.

Python utility functions are provided to analyze these outputs in the pygor.counters submodule.

### Best scenarios
*Output the N best scenarios for each sequence*

Use command `--scenarios N`

The output of this Counter is a semicolon separated values file with one field for each event realization, associated mismatches/errors/mutations indices on the read, the scenario rank, its associated probability and the sequence index. 


### Generation probability
*Estimates the probability of generation of the error free/unmutated ancestor sequence*
By default only outputs an estimator of the probability of generation of the ancestor sequence underlying each sequencing read. See [IGoR's paper][igor_bioarxiv] for details. 

Use command `--Pgen`

### Coverage
*Counts for each genomic nucleotide how many times it has been seen and how many times it was mutated/erroneous*

Use command `--coverage`

## Sequence generation

Using a recombination model and its associated probabilities IGoR can generate random sequences mimicking the raw product of the V(D)J recombination process.

### Sequence generation commands

Reached using the command `-generate N` where *N* is the number of sequences to be generated. The number of sequences to generate must be passed before optional arguments. Optional parameters are the following:

| Command line argument | Description                    |
| :------------- | :------------------------------ |
| `--noerr`  | Generate sequences without sequencing error (the rate and the way those errors are generated is controlled by the model error rate) |
| `--CDR3` | Outputs nucleotide CDR3 from generated sequences. The file contains three fields: CDR3 nucleotide sequence, whether the CDR3 anchors were found (if erroneous/mutated) and whether the sequence is inframe or not. Gene anchors are not yet defined for all the default models shipped with IGoR, use `-set_CDR3_anchors` to set them. |
| `--name myname`  | Prefix for the generated sequences filenames. **Note that setting the *batchname* will change the generated sequences folder name, while setting *--name* will change the file names.** |
| `--seed X`  | Impose *X* as a seed for the random sequence generator. By default a random seed is obtained from the system. |

## Command examples
First as a sanity check try and run the demo code (this will run for a few minutes on all cores available):
```
#!bash
./igor -run_demo
```

Here we give an example with a few commands illustrating a typical workflow. In this example we assume to be executing IGoR from the directory containing the executable.

```
#!bash
WDPATH=/path/to/your/working/directory #Let's define a shorthand for the working directory

#We first read the sequences contained in a text file inside the demo folder
#This will create the align folder in the working directory and the mydemo_indexed_seqs.csv file.
./igor -set_wd $WDPATH -batch foo -read_seqs ../demo/murugan_naive1_noncoding_demo_seqs.txt 

#Now let's align the sequences against the provided human beta chain genomic templates with default parameters
#This will create foo_V_alignments.csv, foo_D_alignments.csv and foo_J_alignments.csv files inside the align folder.
./igor -set_wd $WDPATH -batch foo -species human -chain beta -align --all

#Now use the provided beta chain model to get the 10 best scenarios per sequence
#This will create the foo_output and foo_evaluate and the corresponding files inside
./igor -set_wd $WDPATH -batch foo -species human -chain beta -evaluate -output --scenarios 10

#Now generate 100 synthetic sequences from the provided human beta chain model
#This will create the directory bar_generate with the corresponding files containing the generated sequences and their realizations
./igor -set_wd $WDPATH -batch bar -species human -chain beta -generate 100

```
Since all these commands use several time the same arguments here is some syntactic sugar using more Bash syntax for the exact same workflow with a lighter syntax:
 
```
#!bash
WDPATH=/path/to/your/working/directory #Let's define a shorthand for the working directory
MYCOMMANDS=./igor -set_wd $WDPATH 

$MYCOMMANDS -batch foo -read_seqs ../demo/murugan_naive1_noncoding_demo_seqs.txt #Read seqs
MYCOMMANDS=$MYCOMMANDS -species human -chain beta #Add chain and species commands
$MYCOMMANDS -batch foo -align --all #Align
$MYCOMMANDS -batch foo -evaluate -output --scenarios 10 #Evaluate
$MYCOMMANDS -batch bar -generate 100 #Generate

```

[comment]: # (### Subsampling)

[comment]: # (### Alignments)

[comment]: # (### Inference and Evaluation)

[comment]: # (### Generation)


# Advanced usage
The set of command lines above allows to use predefined models or their topology to study a new dataset. Additionnaly the user can define new models directly using the model parameters file interface. For instance, in order to investigate a conditionnal dependence between two recombination events, the user can simply add or remove an edge in the graph following the syntax defined earlier.

In order to change the set of realizations associated with an event the user can also directly modify a recombination parameters file. Adding or removing realizations should be done with great care as IGoR will use the associated indices to read the corresponding probabilities on the probability array. These indices should be contiguous ranging from 0 to the (total number of realizations -1). Any change in these indices will make the corresponding model marginals file void, and a new one should be automatically created by passing only the model parameters filename to the `-set_custom_model` command. 

Note that changing the GeneChoice realizations can be done automatically (without manually editing the recombination parameter file) by supplying the desired set of genomic templates to IGoR using the `-set_genomic` command. This could be used e.g to define a model for a chain in a species for which IGoR does not supply a model starting from of model for this chain from another species. 

# C++
Although a few command line options are supplied for basic use of IGoR, its full modularity can be used through high level C++ functions on which all previous command lines are built. A section of the main.cpp file is dedicated to accept user supplied code and can be executed using `-run_custom` command line when launching IGoR from the shell. An example of the high level workflow is given in the *run demo* section and the full Doxygen generated documentation is available as PDF. For any question please contact us.

Good practice would be to append the C++ code in the main in the scope where "//Write your custom procedure here" is written. This part of the code is reachable using the `-run_custom` command line argument. This is done so that even after appending some custom code the command line interface is still usable.

# Python
A set of Python codes are shipped with Igor in order to parse IGoR's outputs (alignments,models etc) as the pygor module.

# Contribute

* Your feedbacks are valuable, please send your comments about usability, bug reports and new features you would like to see  
* Code contribution: IGoR was designed to be modular and evolve, please get in touch if you would like to do something new with your data and would like some more guidance on the code structure


# Contact 

For any question please email <quentin.marcou@lpt.ens.fr>
