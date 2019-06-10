## set working directory to the directory containing this script


print("generate synthetic data sets")
source("varselection.generatedata.run.R")

print("generate long MCMC chain")
source("varselection.mcmc.run.R")


## The rest requires running batch scripts.

## The tables were obtained by running varselection.script.R using GNU parallel.
## The script takes 8 arguments:
# JOB_ID
# NRUNS
# design  # 0 for independent design, 1 for correlated
# n       # number of rows
# p       # number of columns
# SNR     # Signal to noise ratio, e.g. 0.5, 1, or 2
# k       # set to zero if only interested in meeting times
# m       # set to zero if only interested in meeting times

# for instance to run the script on 10 machines, 10 times on each machine,
# with either design=0 and design=1, with n=500, with p=1000 and p=5000, and with k=m=0

# parallel -j10 'Rscript varselection.script.R {}' ::: {1..10} ::: 10 ::: {0,1} ::: 500 ::: {1000,5000} ::: {0.5,1,2} ::: 0 ::: 0

## Once this has been run with the desired values of n and p
print("producing tables")
source("varselection.tables.R")

## The figures were obtained as follows.

## For the impact of dimension
print("producing results on the impact of dimension")
source("varselection.differentp.R")
source("varselection.differentp.plots.R")


## For the impact of the hyperparameter kappa
## First, the script "varselection.clusterscript.kappas.R" was run on a cluster.
## This was done by creating a script, say "run.odyssey.varselection.kappas.sh", which contains the following lines

# #!/bin/bash
# #SBATCH -J varselection_rep                     # Job name
# #SBATCH -n 1                                    # Number of cores
# #SBATCH -N 1                                    # All cores on one machine
# #SBATCH -t 0-24:00                              # Runtime in D-HH:MM
# #SBATCH -p shared                               # Partition to submit to (general, shared, serial_requeue)
# #SBATCH --mem-per-cpu=2000M                     # Memory pool for all cores (see also --mem-per-cpu)
# #SBATCH --mail-type=END                         # Type of email notification- BEGIN,END,FAIL,ALL
# #SBATCH --mail-user=username@provider.com       # Email
# #SBATCH --array=1-600                           # Requesting 100 jobs
#
# ## LOAD SOFTWARE ENV ##
# source new-modules.sh
# module load R/3.4.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
# input=varselection.differentkapps.R                                  #****** modify this line
# cd /n/home12/pjacob/                                     #****** modify this line
#
# srun R CMD BATCH --no-save $input out/$input.$SLURM_ARRAY_TASK_ID.out     #****** modify this line


## and then running in a command line
## sbatch run.odyssey.varselection.kappas.sh

## Once these are produced, it remains to run long MCMC chains for comparison,
## and to produce the plots. This is done as follows.
print("producing results on the impact of hyperparameter kappa")
source("varselection.differentkappas.R")
source("varselection.differentkappas.plots.R")
