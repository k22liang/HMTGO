# HMTGO

An R package for a hidden Markov tree model for testing multiple hypotheses corresponding to Gene Ontology gene sets (HMTGO).

### Overview

Testing predefined gene categories has become a common practice for scientists analyzing high throughput transcriptome data.  A systematic way of testing gene categories leads to testing hundreds of null hypotheses that correspond to nodes in a directed acyclic graph.  The relationships among gene categories induce logical restrictions among the corresponding null hypotheses.  An existing fully Bayesian method is powerful but computationally demanding.

We develop a computationally efficient method based on a hidden Markov tree model (HMTM). Our method is several orders of magnitude faster than the existing fully Bayesian method. Through simulation and an expression quantitative trait loci study, we show that the HMTM method provides more powerful results than other existing methods that honor the logical restrictions.

The main function is `HMT`, see the corresponding help file for details.

Citation:
Kun Liang, Chuanlong Du, Hankun You, and Dan Nettleton (2018).  A hidden Markov tree model for testing multiple hypotheses corresponding to gene ontology gene sets. BMC Bioinformatics, 19 107.

### Installation

To install the package in R, use the commands
`````````
library(devtools)
install_github('k22liang/HMTGO', subdir="HMTGO/")
`````````





