## set working directory to the directory containing this script

print("runs on pump failure example")

print("generate estimators for different k and m")
source("pumpfailures.tuning.run.R")

print("generate estimators")
source("pumpfailures.run.R")

print("generate long MCMC run")
source("pumpfailures.mcmc.run.R")

print("regeneration approach")
source("pumpfailures.mykland.run.R")

print("produce plots and tables")
source("pumpfailures.plots.R")
