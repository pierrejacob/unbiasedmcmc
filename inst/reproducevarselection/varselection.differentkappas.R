library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
library(tidyr)
library(viridis)
library(dplyr)
registerDoParallel(cores = detectCores()-2)
#

n <- 500
p <- 1000
SNR <- 1
#
# load data
load(paste0("../varselection.dataSNR", SNR, ".RData"))
# subset data to desired size
Y <- Y[1:n]
X <- X[1:n,1:p]
Y2 <- (t(Y) %*% Y)[1,1]
g <- p^3
#
s0 <- 100
proportion_singleflip <- 0.5
# kappas
nkappas <- 3
kappas <- c(0.1, 1, 2)
# for this kappa
df.mcmc <- data.frame()
# Standard MCMC
nmcmc <- 1e6
burnin <- 1e5
mcmcfilepath <- paste0("varselection.mcmc.differentkappas.n", n, ".p", p, ".RData")
for (ikappa in seq_along(kappas)){
# ikappa <- 1
  print(ikappa)
  kappa <- kappas[ikappa]
  # load model
  vs <- get_variableselection(Y,X,g,kappa,s0,proportion_singleflip)
  prior <- vs$prior
  marginal_likelihood <- vs$marginal_likelihood
  rinit <- vs$rinit
  single_kernel <- vs$single_kernel
  df.mcmc_ <- foreach(irep = 1:10, .combine = rbind) %dorng% {
    current_gamma <- rinit()
    current_pdf <- marginal_likelihood(current_gamma) + prior(current_gamma)
    # chain <- matrix(nrow = nmcmc, ncol = p)
    pdfs <- rep(0, nmcmc)
    sumchain <- rep(0, p)
    # chain[1,] <- current_gamma
    pdfs[1] <- current_pdf
    for (imcmc in 2:nmcmc){
      result <- single_kernel(current_gamma, current_pdf)
      current_gamma <- result$state
      current_pdf <- result$pdf
      # chain[imcmc,] <- current_gamma
      if (imcmc >= burnin){
        sumchain <- sumchain + current_gamma
      }
      pdfs[imcmc] <- current_pdf
    }
    #
    # postburnin <- chain[burnin:nmcmc,]
    postmean <- sumchain/(nmcmc-burnin)
    data.frame(ikappa = ikappa, kappa = kappa, irep = irep, postmean = t(postmean))
  }
  # sapply(1:20, function(i) coda::effectiveSize(chain[burnin:nmcmc,i]))
  # sapply(1:20, function(i) coda::spectrum0(chain[burnin:nmcmc,i])$spec)
  df.mcmc <- rbind(df.mcmc, df.mcmc_)
  save(df.mcmc, file = mcmcfilepath)
}

load(file = mcmcfilepath)

#
getcomponent <- function(x) as.numeric(strsplit(x, "[.]")[[1]][2])
# # getcomponent("postmean.19")
#
dfsub.mcmc <- df.mcmc[,1:23] %>% gather(component, postmean, postmean.1:postmean.20)
dfsub.mcmc$numcomponent <- sapply(dfsub.mcmc$component, getcomponent)
dfsub.mcmc %>% head
#
gmcmc <- ggplot(dfsub.mcmc, aes(x = numcomponent, y = postmean, colour = factor(kappa), group = interaction(kappa, irep))) + geom_point() + geom_line() +
  scale_color_viridis(name = "kappa: ", discrete = TRUE) + xlab("variable") + ylab("inclusion probability")

gmcmc

## load results from the cluster

datafiles <- list.files(path = "output/", pattern = "varselection")
datafiles[1]
df_full <- data.frame()
# compute on the fly quantities required for final estimate
deadline <- 2 * 60 * 60
df_summary <- data.frame()
nkappas <- 3
kappas <- c(0.1, 1, 2)
get_ikappa <- function(igrid) 1 + (igrid-1) %% nkappas

for (ifile in seq_along(datafiles)){
  load(file = paste0("output/varselection.ikappa", get_ikappa(ifile), ".ID", ifile, ".RData"))
  starttime <- c(0, cumsum(durations_))[1:(length(durations_))]
  endtime <- cumsum(durations_)
  df_full <- rbind(df_full, data.frame(rep = rep(ifile, nsamples_),
                                       isample = 1:nsamples_,
                                       durations = durations_,
                                       ikappa = get_ikappa(ifile),
                                       kappa = kappas[get_ikappa(ifile)],
                                       starttime = starttime,
                                       endtime = endtime,
                                       meetings = sapply(results_, function(x) x$meetingtime)))
  # by endtime, we have
  ncompleted <-  which(endtime > deadline)[1] - 1
  if (is.na(ncompleted)){
    ncompleted <- nsamples_
  }
  if (ncompleted == 0){
    ncompleted <- 1
  }
  uestimators <- rep(0, 20)
  uestimators_squared <- rep(0, 20)
  meanduration <- 0
  for (isample in 1:ncompleted){
    uestimators <- uestimators + results_[[isample]]$uestimator[1:20]
    uestimators_squared <- uestimators_squared + (results_[[isample]]$uestimator[1:20])^2
    meanduration <- meanduration + durations_[isample]
  }
  # average estimators produced by this processor
  uestimators <- uestimators / ncompleted
  # average squared estimators
  uestimators_squared <- uestimators_squared / ncompleted
  # average cost
  meanduration <- meanduration / ncompleted
  df_summary <- rbind(df_summary, data.frame(rep = ifile, ikappa = get_ikappa(ifile),
                                             kappa = kappas[get_ikappa(ifile)], ncompleted = ncompleted,
                                             meanduration = meanduration,
                                             uestimators = matrix(uestimators, nrow = 1),
                                             uestimators_squared = matrix(uestimators_squared, nrow = 1)))
}

df_full %>% head

# how long did the first sample take? (in minutes)
summary((df_full %>% filter(isample == 1))$endtime) / 60
# how long has the whole thing been running?
max(df_full$endtime)/3600
# subset df by completion time in seconds
df <- df_full %>% filter(endtime < deadline)
df %>% head
# df[,1:10] %>% head(20)


library(viridis)
g <- ggplot(df, aes(y = rep, yend = rep, x = starttime/60+0.01, xend = endtime/60-0.01, colour = isample)) + geom_segment(lineend = "round", alpha = 1)
g <- g + geom_point() +  scale_color_viridis(name = "sample index:", discrete = FALSE)
g <- g + theme(legend.position = "none") + xlab("time (minutes)") + ylab("processor")  + geom_vline(xintercept = max(df$endtime)/60, linetype = 2)
g

# number of samples produced per processor
summary((df %>% group_by(rep) %>% summarise(max = max(isample)))$max)
# hist((df %>% group_by(rep) %>% summarise(max = max(isample)))$max)
# total number of estimators
nrow(df)

# ggplot(df, aes(x = meetings, y = durations/60, colour = rep)) + geom_point() + ylab("minutes") +
#   scale_color_viridis(discrete = FALSE) + theme(legend.position = "none")

# meeting times
head(df_full)
df_full %>% group_by(ikappa,kappa) %>% summarise(mean(meetings))
max(df_full$meetings)
quantile(df_full$meetings, probs = 0.99)
# meetings <- (df_full %>% filter(isample <= 10) %>% select(meetings))$meetings
qplot(x = df_full$meetings, geom = "blank") + geom_histogram(aes(y = ..density..))
# cat("average meeting time:", mean(df_full$meetings), "\n")
# compute average per processor, for each beta

df_summary %>% head
ncol(df_summary)

# summary((df_summary %>% filter(ikappa == 1))$uestimators.10)
# summary((df_summary %>% filter(ikappa == 1))$uestimators_squared.10)
# df_summary %>% filter(uestimators_squared.10 > 10)
# df_full %>% filter(rep == 445)

df_meanduration <- df_summary %>% group_by(ikappa,kappa) %>% summarise(meanduration = mean(meanduration))
df_meanduration
df_mean <- df_summary %>% group_by(ikappa,kappa) %>% summarise_at(.vars = names(.)[6:25], .funs = c(mean = "mean"))
df_meansquared <- df_summary %>% group_by(ikappa,kappa) %>% summarise_at(.vars = names(.)[26:45], .funs = c(mean = "mean"))
df_sd <- df_summary %>% group_by(ikappa,kappa) %>% summarise_at(.vars = names(.)[6:25], .funs = c(sd = "sd"))

# put all elements in a dataframe for easy plotting
df.plot <- data.frame()
for (ik in 1:3){
  for (ivariable in 1:20){
    estimate <- as.numeric((df_mean %>% filter(ikappa == ik))[paste0('uestimators.', ivariable, '_mean')])
    estimate2 <- as.numeric((df_meansquared %>% filter(ikappa == ik))[paste0('uestimators_squared.', ivariable, '_mean')])
    meanduration <- as.numeric((df_meanduration %>% filter(ikappa == ik))['meanduration'])
    sigma2hat <- meanduration * (estimate2 - estimate^2)
    sderror <- sqrt(sigma2hat) / sqrt(200*deadline) # this one corresponds to Eq (3.5)
    sderror.eq34 <- sqrt(as.numeric((df_sd %>% filter(ikappa == ik))[paste0('uestimators.', ivariable, '_sd')])) / sqrt(200)
    df.plot <- rbind(df.plot, data.frame(ikappa = ik, ivariable = ivariable, kappa = kappas[ik], estimate = estimate, sderror = sderror,
                                         sderror.eq34 = sderror.eq34))
  }
}
df.plot %>% head

gumcmc <- ggplot(df.plot, aes(x = ivariable, y = estimate, colour = factor(kappa), shape = factor(kappa))) + geom_point(size = 2)
gumcmc <- gumcmc + scale_color_viridis(name = "kappa: ", discrete = TRUE) + xlab("variable") + ylab("inclusion probability") + scale_shape(name = "kappa: ")
gumcmc <- gumcmc + geom_errorbar(aes(ymin = estimate - 2 * sderror.eq34, ymax = estimate + 2 * sderror.eq34))
gumcmc


gumcmc2 <- gumcmc + geom_line(data = dfsub.mcmc %>% filter(kappa %in% c(0.1, 1, 2)) %>% select(-component) %>% setNames(c("ikappa", "kappa", "irep", "mcmcestimate", "ivariable")),
                               aes(y = mcmcestimate, group = interaction(kappa, irep)), linetype = 1, alpha = 0.25)
gumcmc2
ggsave(filename = "varselection.estimates.differentkappa.pdf", plot = gumcmc2, width = 8, height = 6)

# gumcmc <- aes(x = numcomponent, y = mean, colour = factor(kappa), group = factor(kappa), shape = factor(kappa)))
# gumcmc <- gumcmc + geom_point(size = 2) + geom_line() +
#   scale_color_viridis(name = "kappa: ", discrete = TRUE) + xlab("variable") + ylab("inclusion probability") + scale_shape(name = "kappa: ")
# gumcmc <- gumcmc + geom_errorbar(aes(ymin = mean - 2 * sd / sqrt(nrep), ymax = mean + 2 * sd / sqrt(nrep)))
# gumcmc <- gumcmc + scale_x_continuous(breaks = c(1,10,20))
# gumcmc
#
# meanduration <- mean(df_summary$meanduration)
# averagesquare <- as.numeric(colMeans(df %>% group_by(rep) %>% summarise_at(.vars = names(.)[7:22], .funs = function(x) mean(x^2)))[2:17])
# sigma2hat <- mean_durations * (averagesquare - estimators^2)
# so compute width of estimators as
# max(df$endtime)
# shat <- sqrt(sigma2hat) / sqrt(500 * max(df$endtime))



