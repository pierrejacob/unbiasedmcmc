## set working directory to the directory containing this script

print("obtain meeting times for single-site Gibbs, on a grid of values of beta")
source("ising.gibbs.meetings.run.R")

print("obtain meeting times for parallel tempering, on a grid of values of beta")
source("ising.swap.meetings.run.R")

print("long run of parallel tempering")
source("ising.mcmc.run.R")

## to run the rest of the calculation you need to run "run.batch.ising.R" many times
## which will populate an output/ subfolder

## See explanations in "run.batch.ising.R" and "run.odyssey.ising.sh".

## Example of command to run the script 10 times, using GNU Parallel and 5 processors:
## parallel -j5 'R CMD BATCH --no-save "--args {1}" run.batch.ising.R out/meeting.job{1}.out'  ::: {1..10}

## Once the batch scripts are run, we can proceed and create the plots.

print("produce plots")
source("ising.plots.R")

