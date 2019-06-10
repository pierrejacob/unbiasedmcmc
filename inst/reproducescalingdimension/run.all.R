
setwd("~/Dropbox/PolyaGamma/code/debiasedmcmc/inst/reproducescalingdimension/")

print("runs on scaling of meeting times with respect to the dimension of a Normal target")

print("produce random walk MH with maximal coupling runs")
source("scaling.rwmh.maximalcoupling.R")
print("produce random walk MH with reflection-maximal coupling runs")
source("scaling.rwmh.reflectionmaxcoupling.R")
print("produce Gibbs sampler runs")
source("scaling.gibbs.R")
print("produce Hamiltonian Monte Carlo runs")
source("scaling.hmc.R")
print("produce plots")
source("scaling.plots.R")
