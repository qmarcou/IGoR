Following are change highlights associated with official releases.  Important
bug fixes are all mentioned, but some internal enhancements are omitted here for
brevity.  Much more detail can be found in the git revision history:

    https://github.com/qmarcou/IGoR

* 1.4.0 (April 11, 2020)

  A brief description of the release content's

  New features:
  - Add feature extracting CDR3s based on alignments (by @alfaceor)
  - New IGL and IGK models for humans (by @alfaceor)

  Bug fixes:
  - fix out of bound access for Dinucleotide markov model internal probabilities (by @Thopic)

  Miscellaneous :
  - New script igor-compute_pgen to calculate pgen for a single sequence (by @alfaceor)
  - Add error messages for unknown supplied chain or output subargument (by @qmarcou)
  - Add new option -getdatadir to output the path of default models directory (by @alfaceor)

* 1.3.0 (August 4, 2018)

  A brief description of the release content's

  New features:
  - Add ---best_gene_only boolean for alignments
  - Allow for template specific alignment offsets
  - Allow for reversed offsets use via command line
  - Pygor package installable via pip (by @penuts7644)
  - New documentation website available at: https://qmarcou.github.io/IGoR/

  Bug fixes:
  - Fix autogen creation of man and html README
  - Fix semi colon separated files for gene anchors (human BCR, mice TRB)
  - Fix pygor package broken imports and functions (by @penuts7644)
  - Fix IGoR directory creation and installation problems (by @smoe)
  - Fix bug on number of scenarios for the scenario counter CL
  - Fix alignment threshold for --ntCDR3 sequences

  Miscellaneous :
  - Split documentation in several files.
  - Create a proper manpage.
  - Improve error handling for alignment issues (report faulty genomic template)
  - Improve error handling for some input files reading
  - make the autogen.sh and build_release scripts executable (by @NickEngland)
  - Clean and make the Pygor package PEP8 compliant (by @penuts7644)

* 1.2.0 (April 18, 2018)

  A brief description of the release content's

  New features:
  - Enable CDR3 sequences alignment using defined gene anchors via --ntCDR3
  - better quality random seed generation
  - show a progress bar for alignment, inference/evaluation and generation.
  - add independent Nmer hypermutation models in the repository
  - add mouse beta chain model

  Bug fixes:
  - replace the use of default_random_engine by a 64bits mersenne twister

  Miscellaneous :
  - make a proper licensing under GNU-GPLv3
  - create a ChangeLog for better traceability
  - clean some Eclipse config file from the repo
  - remove events initialization message duplication
  - transition to asciidoc documentation
  - add juman BCR CDR3 anchors
