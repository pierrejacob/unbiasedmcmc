# set working directory to folder containing this file


print("run various coupled MCMC")
source("germancredit.run.R")

print("run long MCMC")
source("germancredit.mcmc.run.R")

print("produce plots")
source("germancredit.plots.R")

