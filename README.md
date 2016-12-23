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

+ GSL library : currently working on shipping it with IGoR
+ jemalloc (optional although recommended for full parallel proficiency)

# Install

## Linux
Widely tested on several Debian related distros

## MacOS


* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

# Command line tools
Although the full flexibility of IGoR is reachable through C++ highlevel functions (*see next section*) we provide some command line options to perform most frequent tasks on immune receptor sequences.

## Alignments

## Inference

## Outputs 
# C++

# Contribute

* Your feedbacks are valuable, please send your comments about usability and new features you would like to see  
* Code contribution: IGoR was designed to be modular and evolve, please get in touch if you would like to do something new with your data and would like some more guidance on the code structure


# Contact 

For any question please file an issue on github or email quentin.marcou@lpt.ens.fr