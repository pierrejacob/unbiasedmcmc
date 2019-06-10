# debiasedmcmc

This package contains scripts to reproduce the results of the arXiv report "Unbiased Markov chain Monte Carlo with couplings", by Pierre E. Jacob, John O'Leary and Yves F Atchade.

This is not a general-purpose statistical software, but a collection of scripts intended to reproduce figures and tables. Use at your own risk!

The Polya-Gamma sampler is taken from the package BayesLogit of Nick Polson, James Scott, and Jesse Windle.
https://cran.r-project.org/src/contrib/Archive/BayesLogit/
That package can be downloaded and installed with "R CMD INSTALL".
Or in R via
packageurl <- "https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
Note that the above requires gfortran. For instance on Mac OS X this can be found on https://github.com/fxcoudert/gfortran-for-macOS/releases.

The folder inst/reproduce/ contains the scripts to reproduce the figures. Each sub-folder has a "run all" script to run all the scripts in the correct order.
The folder inst/check/ contains internal checks, and vignettes/ contains tutorials on
how to use the package's main functions.
