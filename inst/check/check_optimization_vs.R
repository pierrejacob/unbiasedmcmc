library(debiasedmcmc)
setwd("~/Dropbox/UnbiasedMCMCResults/reproduce/")
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores() - 1)
#

# simulate data
n <- 500
p <- 10000
SNR <- 3
s_star <- 10
s0 <- 100
sigma0 <- 1
beta_star <- SNR * sqrt(sigma0^2 * log(p) / n) * c(2,-3,2,2,-3,3,-2,3,-2,3, rep(0, p-10))
# independent design
X <- matrix(rnorm(n * p), nrow = n, ncol = p) # fast_rmvnorm_chol(n, rep(0, p), diag(1, p, p))
X <- scale(X)
Y <- X %*% matrix(beta_star, ncol = 1) + rnorm(n, 0, sigma0)
Y <- scale(Y)
Y2 <- (t(Y) %*% Y)[1,1]
g <- p^3
kappa <- 1
proportion_singleflip <- 0.5
vs <- get_variableselection(Y,X,g,kappa,s0,proportion_singleflip)
prior <- vs$prior
marginal_likelihood <- vs$marginal_likelihood
rinit <- vs$rinit
single_kernel <- vs$single_kernel
coupled_kernel <- vs$coupled_kernel

nmcmc <- 2000
gammas <- matrix(nrow = nmcmc, ncol = p)
pi_gammas <- rep(0, nmcmc)
gamma <- rinit()
pi_gamma <- vs$prior(gamma) + vs$marginal_likelihood(gamma)


cppFunction('
double marginal_likelihood_test_(Eigen::VectorXf selection, const Eigen::MatrixXf & X, const Eigen::VectorXf & Y, double Y2, double g){
// std::cerr << EIGEN_WORLD_VERSION << ". " << EIGEN_MAJOR_VERSION << ". " << EIGEN_MINOR_VERSION << std::endl;
  double l = 0.;
  int n = X.rows();
  int p = X.cols();
  int s = selection.sum();
  if (s > 0){
    Eigen::MatrixXf Xselected(n,s);
    int counter = 0;
    for (int column = 0; column < p; column++){
      if (selection(column)){
        Xselected.col(counter) = X.col(column);
        counter++;
      }
    }
    Eigen::BDCSVD<Eigen::MatrixXf > svd(Xselected, Eigen::ComputeThinU);
    Eigen::VectorXf z = svd.matrixU().transpose() * Y.head(s);
    l = z.dot(z);
//    l = Y.transpose() * Xselected * ((Xselected.transpose() * Xselected).inverse()) * Xselected.transpose() * Y;
  } else {
    l = 0.;
  }
  l = -((double) s + 1.) / 2. * log(g + 1.) - (double) n / 2. * log(Y2 - g / (g + 1.) * l);
  return l;
}', depends = "RcppEigen")

vs$marginal_likelihood(gamma)
marginal_likelihood_test_(gamma, X, Y, Y2, g)

library(microbenchmark)
microbenchmark(vs$marginal_likelihood(rinit()),
               marginal_likelihood_test_(rinit(), X, Y, Y2, g))

#
#
#
#
# for (imcmc in 1:nmcmc){
#   result <- single_kernel(gamma, pi_gamma)
#   gamma <- result$state
#   pi_gamma <- result$pdf
#   gammas[imcmc,] <- gamma
#   pi_gammas[imcmc] <- pi_gamma
# }
#
# plot(pi_gammas)
#
# posterior_mean <- colMeans(gammas[10000:nmcmc,])
# plot(posterior_mean)
# posterior_mean[1:20]
# sum(posterior_mean > 1e-5)

