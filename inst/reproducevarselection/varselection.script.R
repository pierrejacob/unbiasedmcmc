# Script that sets an initial seed using rlecuyer's package
# and compute a number of estimators using that seed
#
arguments <- commandArgs(TRUE)
if (length(arguments)!=7){
  cat("needs 7 integers as arguments; JOB_ID, NRUNS, n, p, SNR, k, m")
  q()
}

# arguments <- c("1", "10", "1000", "500", "3", "0", "0")
#
JOB_ID <- as.numeric(arguments[1])
NRUNS <- as.numeric(arguments[2]) # number of runs for this job
n <- as.numeric(arguments[3])
p <- as.numeric(arguments[4])
SNR <- as.numeric(arguments[5])
k <- as.numeric(arguments[6])
m <- as.numeric(arguments[7])

# file paths
workingdirectory <- "~/Dropbox/UnbiasedMCMCResults/varselection/"
resultpath <- paste0("varselection.SNR", SNR, ".n", n, ".p", p, ".k", k, ".m", m, ".job", JOB_ID, ".RData")
#
setwd(workingdirectory)
# load packages
library(rlecuyer)
library(debiasedmcmc)
setmytheme()
# initial seed
.lec.SetPackageSeed(c(42, 66, 101, 123454, 7, 54321))
nstream <- 1000 # number larger than total # of processors expected to run this script
stream.names <- paste(1:nstream)
.lec.CreateStream(stream.names)
.lec.CurrentStream(paste(JOB_ID))

### beginning of job
load(paste0("~/Dropbox/UnbiasedMCMCResults/varselection/varselection.dataSNR", SNR, ".RData"))
Ysub <- Y[1:n]
Xsub <- X[1:n,1:p]
Y2 <- (t(Ysub) %*% Ysub)[1,1]
g <- p^3

s0 <- 100
kappa <- 0.1
proportion_singleflip <- 0.5

vs <- get_variableselection(Ysub,Xsub,g,kappa,s0,proportion_singleflip)

prior <- vs$prior
marginal_likelihood <- vs$marginal_likelihood
rinit <- vs$rinit
single_kernel <- vs$single_kernel
coupled_kernel <- vs$coupled_kernel
coupled_chains <- vs$coupled_chains
unbiasedestimator <- vs$unbiasedestimator

result <- list()
for (irun in 1:NRUNS){
  cat("Run #", irun, "\n")
  ue <- unbiasedestimator(single_kernel, coupled_kernel, rinit, h = function(x) x, k = k, m = m)
  result[[irun]] <- list(irun = irun, JOB_ID = JOB_ID,
                         mcmcestimator = ue$mcmcestimator, correction = ue$correction,
                         uestimator = ue$uestimator, meetingtime = ue$meetingtime, niterations = ue$iteration)
  # cc <- coupled_chains(single_kernel, coupled_kernel, rinit, m = m)
  # measure <- debiasedmcmc:::get_measure_(cc, k, m)
  # approximation <- data.frame(rep = rep(irun, length(measure$weights)), weights = measure$weights, atoms = measure$atoms)
  # approximation <- data.frame(debiasedmcmc:::prune_(as.matrix(approximation[do.call(order, as.list(approximation[,3:ncol(approximation)])),])))
  # approximation <- approximation[approximation$weights != 0,]
  # result[[irun]] <- list(irun = irun, JOB_ID = JOB_ID, approximation = approximation, meetingtime = cc$meetingtime, niterations = cc$iteration)
  save(result, file = resultpath)
}
.lec.CurrentStreamEnd()
