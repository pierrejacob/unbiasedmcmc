## set working directory to the directory containing this script

print("run long MCMC")
source("baseball.mcmc.run.R")

print("run coupled chains with m=100")
source("baseball.run.R")

print("produce plots")
source("baseball.plots.R")
