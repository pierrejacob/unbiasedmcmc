#!/bin/bash
#SBATCH -J ising_rep                           # Job name
#SBATCH -n 1                                    # Number of cores
#SBATCH -N 1                                    # All cores on one machine
#SBATCH -t 0-24:00                              # Runtime in D-HH:MM
#SBATCH -p shared                               # Partition to submit to (general, shared, serial_requeue) *************** modify this line
#SBATCH --mem-per-cpu=2000M                     # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END                         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=username@fas.harvard.edu    # Email *************** modify this line
#SBATCH --array=1-100                           # *************** modify this line

## LOAD SOFTWARE ENV ##
source new-modules.sh
module load R/3.4.2-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
input=run.batch.ising.R           # *************** modify this line
cd /n/home12/pjacob/isingmodel/   # *************** modify this line

srun R CMD BATCH --no-save "--args $SLURM_ARRAY_TASK_ID" $input out/$input.$SLURM_ARRAY_TASK_ID.out     #****** modify this line


### How to run the calculations on the Harvard Odyssey cluster?
### command to ssh to the cluster, e.g.
# ssh username@odyssey.rc.fas.harvard.edu

###
# First, install debiasedmcmc package on the cluster (e.g. upload zip file and R CMD INSTALL).
# Do this in an interactive session, e.g.
####
# srun --pty -p test -t 20 --mem 3000 /bin/bash
### to have access to R:
# source new-modules.sh
# module load R/3.4.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
# R CMD INSTALL package.tar.gz


### commands such as the following may be useful when it comes to transferring files
# scp /home/bla*bla blabla@odyssey.rc.fas.harvard.edu:~/blabla/
# scp blabla@odyssey.rc.fas.harvard.edu:~/blabla/output/*RData /home/blabla/
### see more details on: https://www.rc.fas.harvard.edu/resources/documentation/copying-data-to-and-from-odyssey-using-scp/

### Create appropriate folders (subfolders named "output" and "out" in particular).
### tune run.batch.ising.R appropriately (it should save output in output/ subfolder).

### Make sure path names are OK in run.odyssey.ising.sh.
### choose partition, size of array (e.g. 1-500 for 500 jobs), runtime
### To start the jobs:
# sbatch run.odyssey.ising.sh.

### To monitor jobs,
# sacct -j JOBID
### or
# squeue -j 9999999
### Or go to the website https://portal.rc.fas.harvard.edu/jobs/
### To cancel a job,
# scancel 9999999
