IGoR: Inference and Generation Of Repertoires
=============================================


This repository contains all sources and models useful to infer V(D)J
recombination related processes for TCR or BCR sequencing data using
*IGoR*.

IGoR's exhaustive *_documentation can be found http://qmarcou.github.io/IGoR[here]_*!

Download the latest release https://github.com/qmarcou/IGoR/releases[here].

[[quick-summary]]
Quick summary
-------------

IGoR is a C++ software designed to infer V(D)J recombination related
processes from sequencing data such as:

* Recombination model probability distribution
* Hypermutation model
* Best candidates recombination scenarios
* Generation probabilities of sequences (even hypermutated)

The following paper describes the methodology, performance tests and
some new biological results obtained with IGoR:

https://www.nature.com/articles/s41467-018-02832-w[_High-throughput
immune repertoire analysis with IGoR_], _Nature Communications_, (2018)
Quentin Marcou, Thierry Mora, Aleksandra M. Walczak

Its heavily object oriented and modular style was designed to ensure
long term support and evolvability for new tasks in assessing TCR and
BCR receptors features using modern parallel architectures.

IGoR is a free (as in freedom) software released under the
https://www.gnu.org/licenses/quick-guide-gplv3.html[GNU-GPLv3] license.

[[version]]
Version
-------

include::./docs/asciidoc/version.adoc[]
Download the latest release https://github.com/qmarcou/IGoR/releases[here].

Documentation
-------------
A comprehensive documentation of IGoR is available on the internet at: http://qmarcou.github.io/IGoR[http://qmarcou.github.io/IGoR] or in your local installation of IGoR as the `./docs/index.html` file. 

include::./docs/asciidoc/general/contribute.adoc[]


[[contact]]
Contact
-------

For any question or issue please open an
https://github.com/qmarcou/IGoR/issues[issue] or email mailto:quentin.marcou@lpt.ens.fr[us].


Copying
-------
Free use of IGoR is granted under the terms of the https://www.gnu.org/licenses/quick-guide-gplv3.html[GNU General Public License version 3]
(GPLv3).
