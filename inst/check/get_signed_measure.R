library(debiasedmcmc)
library(ggthemes)
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores())

# this script plays with the computation of estimators
# using signed measure representations
logtarget <- function(x) fast_dmvnorm(matrix(x, nrow = 1), mean = c(0,0.5), covariance = diag(c(0.5, 0.2)))
Sigma_proposal <- diag(c(1, 1))
rinit <- function() fast_rmvnorm(1, c(10, 5), diag(c(3,3)))
kernels <- get_mh_kernels(logtarget, Sigma_proposal, dimension = 2)

# logtarget <- function(x) fast_dmvnorm(matrix(x, nrow = 1), mean = c(0), covariance = diag(c(0.5), 1, 1))
# Sigma_proposal <- diag(c(1))
# rinit <- function() fast_rmvnorm(1, c(5), diag(c(3), 1, 1))
# kernels <- get_mh_kernels(logtarget, Sigma_proposal, dimension = 1)

xtest <- rinit()
xtest <- kernels$kernel(xtest)
xtest2 <- rinit()

kernels$coupled_kernel(xtest, xtest2)


nsamples <- 1000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(kernels$kernel, kernels$coupled_kernel, rinit = rinit)
}
sapply(c_chains_, function(x) x$meetingtime) %>% hist

cx <- coupled_chains(kernels$kernel, kernels$coupled_kernel, rinit = rinit)
cx$samples1 %>% head

measure <- debiasedmcmc:::get_measure_(cx, 0, 0)
sum(measure$weights * measure$atoms[,1])
sum(measure$weights * measure$atoms[,2])

H_bar(cx, h = function(x) x, 0, 0)

## now with k,m larger
k <- floor(as.numeric(quantile(sapply(c_chains_, function(x) x$meetingtime), probs = 0.95)))
m <- 5 * k
c_chains_ <- foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(kernels$kernel, kernels$coupled_kernel, rinit = rinit, m = m)
}
approximation <- c_chains_to_dataframe(c_chains_, k, m)

# sum(approximation$weights * approximation$atoms)
# mean(sapply(X = c_chains_, FUN = function(x) H_bar(x, h = function(v) v[1], k, m)))

sum(approximation$weights * approximation$atoms.1)
mean(sapply(X = c_chains_, FUN = function(x) H_bar(x, h = function(v) v[1], k, m)))

sum(approximation$weights * approximation$atoms.2)
mean(sapply(X = c_chains_, FUN = function(x) H_bar(x, h = function(v) v[2], k, m)))

approximation %>% group_by(rep) %>% summarise(sumw = sum(weights)) %>% summary
sum(approximation$rep %>% unique)
nsamples * (nsamples + 1) / 2
#
# approximation <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
#   measure <- debiasedmcmc:::get_measure_(c_chains_[[irep]], k, m)
#   data.frame(rep = rep(irep, length(measure$weights)), weights = measure$weights, atoms = measure$atoms)
# }
# approximation$weights <- approximation$weights / nsamples
# approximation %>% tail
#
# mean(approximation$weights < 0)
# sum(approximation$weights)
#
# # note: we can remove duplicates if we want a smaller data.frame
# cppFunction('DataFrame prune_(const DataFrame & df){
#   NumericVector atoms1 = df["atoms.1"];
#   NumericVector atoms2 = df["atoms.2"];
#   NumericVector weights = df["weights"];
#   NumericVector rep = df["rep"];
#   for (int i=1; i < df.rows(); i++){
#     if (rep(i) == rep(i-1)){
#       if (std::abs(atoms1(i) - atoms1(i-1)) < 1e-30 && std::abs(atoms2(i) - atoms2(i-1)) < 1e-30){
#         weights(i)  += weights(i-1);
#         weights(i-1) = 0;
#       }
#     }
#   }
#   return DataFrame::create(Named("rep") = rep, Named("atoms.1") = atoms1,
#                            Named("atoms.2") = atoms2, Named("weights") = weights);
#   }')
#
# cppFunction('NumericMatrix prune_2(const NumericMatrix & df){
# NumericMatrix copieddf = clone(df);
# double s = 0;
# for (int i=1; i < copieddf.rows(); i++){
#   if (copieddf(i,0) == copieddf(i-1,0)){
#     s = 0;
#     for (int j = 2; j < copieddf.cols(); j++){
#       s += std::abs(copieddf(i,j) - copieddf(i-1,j));
#     }
#     if (s < 1e-20){
#       copieddf(i, 1)  += copieddf(i-1, 1);
#       copieddf(i-1, 1) = 0.;
#     }
#   }
# }
# return copieddf;
#             }')
#
# # test_(as.matrix(approximation %>% arrange(atoms.1)))
# approximation.df2 <- data.frame(prune_2(as.matrix(approximation  %>% arrange(atoms.1))))
# approximation.df2 <- approximation.df2[approximation.df2$weights != 0,]
# approximation.df2 %>% head
#
# prune <- function(approximation.df){
#   sorted.df <- approximation.df %>% arrange(atoms.1, atoms.2)
#   sorted.df <- prune_(sorted.df)
#   return(sorted.df[sorted.df$weights != 0,])
# }
#
# approximation.pruned <- prune(approximation)
#
# sum(approximation.pruned$weights * approximation.pruned$atoms.2)
# sum(approximation.df2$weights * approximation.df2$atoms.2)
# mean(sapply(X = c_chains_, FUN = function(x) H_bar(x, h = function(v) v[2], k, m)))

