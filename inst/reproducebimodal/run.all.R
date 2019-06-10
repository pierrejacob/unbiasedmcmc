## set working directory to the directory containing this script

print("runs on bimodal target")

print("generate meeting times for different MCMC proposals / different initial distribution")
source("bimodal.run.R")

print("generate long MCMC run")
source("bimodal.mcmc.run.R")

print("produce plots and tables")
source("bimodal.plots.R")
