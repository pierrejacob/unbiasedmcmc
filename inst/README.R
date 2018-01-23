### This file explains how to reproduce the figures of the article

# The files of the debiasedmcmc/inst/reproduce*
# contain the scripts to reproduce the figures, model by model.
#
# The files finishing in .run.R are the ones creating the result files;
# these scripts take several hours to run on a desktop computer. One should specify a
# working directory before running the files; the results are saved in the current working directory.
#
# The files finishing in .plots.R create the plots, assuming that the .run.R files
# have been run already. Indeed the plot files require specific .RData files that are
# created by the .run.R files.
#
# For the variable selection example, a script named varselection.script.R is in the folder.
# This script is meant to be run via GNU Parallel in the command line; examples
# of commands can be found in varselection.script.examples.run.R
# Note that there is a hardcoded path to change in varselection.script.R !
