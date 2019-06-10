library(debiasedmcmc)
library(dplyr)
library(doParallel)
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores() - 1)
#

## meeting times
SNR <- 1
n <- 500
get_meeting_times <- function(SNR, n, p, k){
  resultsfiles <- list.files(pattern = paste0("varselection.SNR", SNR, ".n", n, ".p", p, ".k", k, ".*"))
  results <- list()
  iresult <- 1
  meetingtimes <- c()
  resultsfiles
  nfiles <- length(resultsfiles)
  for (ifile in 1:nfiles){
    load(resultsfiles[ifile])
    for (ires in 1:length(result)){
      results[[iresult]] <- result[[ires]]
      iresult <- iresult + 1
    }
  }
  return(results)
}

results <- get_meeting_times(SNR, n, 1000, 0)
meetings <- sapply(results, function(x) x$meetingtime)
meetings %>% length
results[[1]]$meetingtime

# ps <- c(100, 500, 1000)
ps <- c(100, 250, 500, 750, 1000)
meeting.df <- data.frame()
for (p in ps){
  meetings <- get_meeting_times(SNR, n, p, 0)
  meetings <- sapply(meetings, function(x) x$meetingtime)
  nmeetings <- length(meetings)
  meeting.df <- rbind(meeting.df, data.frame(p = rep(p, nmeetings), meetingtime = meetings))
}


meeting.df %>% group_by(p) %>% summarize(mean = mean(meetingtime), q90 = floor(as.numeric(quantile(meetingtime, probs = 0.9))),
                                         q95 = floor(as.numeric(quantile(meetingtime, probs = 0.95))),
                                         q99 = floor(as.numeric(quantile(meetingtime, probs = 0.99))),
                                         n = n())

g <- ggplot(meeting.df, aes(x = p, y = meetingtime / p)) + geom_point() + ylab("meeting times / p")
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.1)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.2)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.3)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.4)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.5)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.6)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.7)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.8)), geom = "line", linetype = 2)
g <- g + stat_summary(fun.y = function(v) as.numeric(quantile(v, probs = 0.9)), geom = "line", linetype = 2)
g

ggsave(filename = "varselection.meetingtimes.pdf", plot = g, width = 6, height = 5)
## estimation for p = 1000

SNR <- 1
n <- 500
p <- 1000
# k <- 40000
# m <- 200000
k <- 50000
m <- 100000

results <- list()
iresult <- 1
resultsfiles <- list.files(pattern = paste0("varselection.SNR", SNR, ".n", n, ".p", p, ".k", k, "*"))
meetingtimes <- c()
resultsfiles
nfiles <- length(resultsfiles)
for (ifile in 1:nfiles){
  load(resultsfiles[ifile])
  for (ires in 1:length(result)){
    results[[iresult]] <- result[[ires]]
    iresult <- iresult + 1
  }
}
nrep <- results %>% length
nrep

results[[1]]$meetingtime
sapply(results, function(x) x$meetingtime) %>% hist
sapply(results, function(x) x$meetingtime) %>% max

estimators <- foreach(iresult = 1:length(results), .combine = rbind) %dopar% {
  results[[iresult]]$uestimator
}

postmeans <- colMeans(estimators)
sdpostmeans <- apply(estimators, 2, sd)
plot(postmeans[1:20], type = "l", lty = 3, ylim = c(0,1))
lines(postmeans[1:20] + 2 * sdpostmeans[1:20]/sqrt(nrep))
lines(postmeans[1:20] - 2 * sdpostmeans[1:20]/sqrt(nrep))

plot(postmeans[21:p], type = "l", lty = 3, ylim = c(0,1))
lines(postmeans[21:p] + 2 * sdpostmeans[21:p]/sqrt(nrep))
lines(postmeans[21:p] - 2 * sdpostmeans[21:p]/sqrt(nrep))


plot(sdpostmeans, type = "l")
load("varselection.SNR1.n500.p1000.mcmc.RData")
postmeans[1:20]
postmean[1:20]
mcmcvar

df <- data.frame(p = 1:20, mcmcpostmean = postmean[1:20],
           confmin = postmeans[1:20] - 2 * sdpostmeans[1:20]/sqrt(nrep),
           confmax = postmeans[1:20] + 2 * sdpostmeans[1:20]/sqrt(nrep))

g <- ggplot(df, aes(x = p, ymin = confmin, ymax = confmax, xmin = p, xmax = p)) + geom_errorbar(alpha = 1)
g <- g + geom_line(aes(y = mcmcpostmean), linetype = 3)
g <- g + xlab("variable") + ylab("inclusion probability")
g

ggsave(filename = "varselection.estimates.pdf", plot = g, width = 6, height = 5)

