# debiasedmcmc

This package contains scripts to reproduce the results of the arXiv report "Unbiased Markov chain Monte Carlo with couplings", by Pierre E. Jacob, John O'Leary and Yves F Atchade.

The Polya-Gamma sampler is taken from the package BayesLogit of Nick Polson, James Scott, and Jesse Windle.
https://cran.r-project.org/src/contrib/Archive/BayesLogit/
That package can be downloaded and installed with "R CMD INSTALL".
Or in R via
packageurl <- "https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
Note that the above requires libgfortran.

The folder inst/reproduce/ contains the scripts to reproduce the figures. Each sub-folder has a "run all" script to run all the scripts in the correct order.
The folder inst/check/ contains internal checks, and vignettes/ contains tutorials on
how to use the package's main functions.
