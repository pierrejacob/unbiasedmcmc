# load packages
library(unbiasedmcmc)
library(dplyr)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())

#
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
curve(sapply(x, function(v) exp(target(v))), from = -10, to = 10)

load(file = "bimodal.c_chains.easy.RData")

x <- as.numeric(names(table(meetingtimes.easy)))
y <- as.numeric(table(meetingtimes.easy)) / length(meetingtimes.easy)
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g

print(mean(meetingtimes.easy))
print(quantile(meetingtimes.easy, probs = c(0.9, 0.99)))


nsamples <- length(c_chains.easy)
taus <- sapply(c_chains.easy, FUN = function(x) x$meetingtime)
mmax <- c_chains.easy[[1]]$iteration
testfunction <- function(x) x > 3
summary(taus)
# ks <- c(1, 50, 100)
# mfactor <- c(1, 10, 100)
ks <- c(1, 100, 200)
mfactor <- c(1, 10, 20)


ineff.df <- data.frame()
for (ik in 1:length(ks)){
  for (im in 1:length(mfactor)){
    k <- ks[ik]
    m <- mfactor[im] * k
    estimators <-  foreach(irep = 1:nsamples) %dorng% {
      H_bar(c_chains.easy[[irep]], h = testfunction, k = k, m = m)
    }
    meanh <- mean(unlist(estimators))
    v <- var(unlist(estimators))
    c <- mean(2 * taus + pmax(1, m + 1 - taus))
    ineff.df <- rbind(ineff.df, data.frame(k = k, mfactor = paste0(mfactor[im], " times k"), c = c, meanh = meanh, v = v))
  }
}

ineff.df %>% head
ineff.df$inef <- (ineff.df$c * ineff.df$v)
ineff.df

exact_x2 <- integrate(f = function(x) sapply(x, function(v) (v>3) * exp(target(v))), lower = -20, upper = 20, subdivisions = 1e5)
exact_x2$value
load("bimodal.mcmc.RData")
mcmcvar.easy
ineff.df$inef <- ineff.df$inef / mcmcvar.easy

library(xtable)
formatted.df <- ineff.df %>% select(-meanh)
formatted.df$k <- paste0(formatted.df$k)
cap <- "Cost, variance and inefficiency divided by MCMC asymptotic variance, for various choices of $k$ and $m$, for the test function $h:x\\mapsto \\mathds{1}(x > 3)$, in the mixture target example of Section \\ref{subsec:Random-Walk-Metropolis-bimodal}. \\label{table:mixture:easy}"
colnames(formatted.df) <- c("k", "m", "expected cost", "variance", "inefficiency / MCMC")
formatted.df$`expected cost` <- format(formatted.df$`expected cost`, digits = 2)
formatted.df$`variance` <- format(formatted.df$`variance`, digits = 2, scientific = TRUE)
formatted.df$`inefficiency / MCMC` <- format(formatted.df$`inefficiency / MCMC`, digits = 2)
formatted.df
formatted.df <- xtable(formatted.df, caption = cap)
print.xtable(formatted.df, include.rownames = FALSE, include.colnames = TRUE, file = "bimodal.inefficiency.tex")

# ggsave(filename = "bimodal.meetingtime.pdf", plot = g, width = 7, height = 7)
##
nclass <- 50
k <- 200
m <- 2000
histogram <- histogram_c_chains(c_chains.easy, 1, k, m, nclass = nclass)
g <- plot_histogram(histogram, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "bimodal.histogram.easy.pdf", plot = g, width = 5, height = 5)


load(file = "bimodal.c_chains.intermediate.RData")

print(mean(meetingtimes.intermediate))
print(quantile(meetingtimes.intermediate, probs = c(0.9, 0.99)))

x <- as.numeric(names(table(meetingtimes.intermediate)))
y <- as.numeric(table(meetingtimes.intermediate)) / length(meetingtimes.intermediate)
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g
nclass <- 50
k <- 20000
m <- 30000
histogram <- histogram_c_chains(c_chains.intermediate, 1, k, m, nclass = nclass)
g <- plot_histogram(histogram, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "bimodal.histogram.intermediate.pdf", plot = g, width = 5, height = 5)

estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains.intermediate[[irep]], h = testfunction, k = k, m = m)
}
meanh <- mean(unlist(estimators))
v <- var(unlist(estimators))
s <- 2*sqrt(v / length(c_chains.intermediate))
cat(meanh-s, meanh+s)


load(file = "bimodal.c_chains.hard.RData")

print(mean(meetingtimes.hard))
print(quantile(meetingtimes.hard, probs = c(0.9, 0.99)))
max(meetingtimes.hard)

x <- as.numeric(names(table(meetingtimes.hard)))
y <- as.numeric(table(meetingtimes.hard)) / length(meetingtimes.hard)

g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g
nclass <- 50
k <- 50
m <- 500
histogram <- histogram_c_chains(c_chains.hard, 1, k, m, nclass = nclass)
g <- plot_histogram(histogram, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "bimodal.histogram.hard.pdf", plot = g, width = 5, height = 5)

estimators <-  foreach(irep = 1:length(c_chains.hard)) %dorng% {
  H_bar(c_chains.hard[[irep]], h = testfunction, k = k, m = m)
}
meanh <- mean(unlist(estimators))
v <- var(unlist(estimators))
s <- 2*sqrt(v / length(c_chains.hard))
cat(meanh-s, meanh+s)

estimators.extra <-  foreach(irep = 1:length(c_chains.hard.extra)) %dorng% {
  H_bar(c_chains.hard.extra[[irep]], h = testfunction, k = k, m = m)
}

estimators.all <- c(unlist(estimators), unlist(estimators.extra))
meanh <- mean(estimators.all)
v <- var(estimators.all)
s <- 2*sqrt(v / length(estimators.all))
cat(meanh-s, meanh+s)
