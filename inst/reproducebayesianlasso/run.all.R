# set working directory to folder containing this file

print("run various coupled MCMC")
source("bayesianlasso.run.R")

print("run long MCMC")
source("bayesianlasso.mcmc.run.R")
# this is used to get effective sample size as a function of lambda

print("produce plots")
source("bayesianlasso.plots.R")

