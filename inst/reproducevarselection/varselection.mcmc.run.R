library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(1)

# file paths
n <- 500
SNR <- 1
p <- 1000
#
load(paste0("varselection.dataSNR", SNR, ".RData"))

Y <- Y[1:n]
X <- X[1:n,1:p]
Y2 <- (t(Y) %*% Y)[1,1]
g <- p^3

s0 <- 100
kappa <- .1
proportion_singleflip <- 0.5

vs <- get_variableselection(Y,X,g,kappa,s0,proportion_singleflip)

prior <- vs$prior
marginal_likelihood <- vs$marginal_likelihood
rinit <- vs$rinit
single_kernel <- vs$single_kernel

nmcmc <- 2e4

current_gamma <- rinit()
current_pdf <- marginal_likelihood(current_gamma) + prior(current_gamma)
chain <- matrix(nrow = nmcmc, ncol = p)
pdfs <- rep(0, nmcmc)
chain[1,] <- current_gamma
pdfs[1] <- current_pdf
for (imcmc in 2:nmcmc){
  result <- single_kernel(current_gamma, current_pdf)
  current_gamma <- result$state
  current_pdf <- result$pdf
  chain[imcmc,] <- current_gamma
  pdfs[imcmc] <- current_pdf
}
#

plot(pdfs, type = "l")

burnin <- nmcmc/2

postburnin <- chain[burnin:nmcmc,]
postmean <- colMeans(postburnin)

# plot(postmean, ylim = c(0,1), type = "l")
plot(postmean[1:20], ylim = c(0,1), type = "b")
print(postmean[1:20])

library(coda)
mcmcvar <- spectrum0(postburnin[,1:10])
mcmcvar
# save(nmcmc, burnin, pdfs, postmean, mcmcvar, file = "varselection.SNR1.n500.p1000.mcmc.RData")
#
