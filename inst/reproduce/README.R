### This file explains how to reproduce the figures of the article
# The files of the debiasedmcmc/inst/reproduce
# contain the scripts to reproduce the figures.
#
# The files finishing in .run.R are the ones creating the result files;
# these scripts take several hours to run on a desktop computer.
#
# The files finishing in .plots.R create the plots, assuming that the .run.R files
# have been run already. Indeed the plot files require specific .RData files that are
# created by the .run.R files.

# folders
scriptfolder <- "~/path to debiasedmcmc/inst/reproduce"
resultsfolder <- "~/path to a folder where large RData files will be created"
#
setwd(resultsfolder)

### to create the diagram (Figure 1)
source(file.path(scriptfolder, "diagram.R"))

### Baseball batting data example
# run
source(file.path(scriptfolder, "baseball.run.R"))
source(file.path(scriptfolder, "baseball.tuning.run.R"))
source(file.path(scriptfolder, "baseball.mcmc.run.R"))
# plots
source(file.path(scriptfolder, "baseball.plots.R"))

### Pump failure data example
# run
source(file.path(scriptfolder, "pumpfailures.run.R"))
source(file.path(scriptfolder, "pumpfailures.tuning.run.R"))
source(file.path(scriptfolder, "pumpfailures.mcmc.run.R"))
# plots
source(file.path(scriptfolder, "pumpfailures.plots.R"))

### Scaling with the dimension
# run
source(file.path(scriptfolder, "scaling.rwmh.run.R"))
source(file.path(scriptfolder, "scaling.gibbs.run.R"))
# plots
source(file.path(scriptfolder, "scaling.plots.R"))

### Bimodal target example
# run
source(file.path(scriptfolder, "mixture.run.R"))
# plots
source(file.path(scriptfolder, "mixture.plots.R"))

### German credit with PGG (to clean, see files in ../needscleaning)
# run
#...
# plots
#...

### Diabetes with Bayesian Lasso
# run
source(file.path(scriptfolder, "diabetes.p10.run.R"))
source(file.path(scriptfolder, "diabetes.p64.run.R"))
source(file.path(scriptfolder, "diabetes.p64.refine.run.R"))

# plots
source(file.path(scriptfolder, "diabetes.plots.R"))

### Cut distribution on Plummer's example
# run
source(file.path(scriptfolder, "plummer.cut.run.R"))
# plots
source(file.path(scriptfolder, "plummer.cut.plots.R"))


### Cut distribution on Zigler's propensity score example
# run
source(file.path(scriptfolder, "propensity.stage1.run.R"))
source(file.path(scriptfolder, "propensity.stage2.run.R"))
source(file.path(scriptfolder, "propensity.doublemcmc.run.R"))
# plots
source(file.path(scriptfolder, "propensity.plots.R"))


