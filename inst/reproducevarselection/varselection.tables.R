library(unbiasedmcmc)
# setmytheme()
library(dplyr)
library(tidyr)
rm(list = ls())
set.seed(1)

# the function gets the results from all the jobs matching the
# given parameters
get_results <- function(design, SNR, n, p, k){
  resultsfiles <- list.files(pattern = paste0("varselection.design", design, ".SNR", SNR, ".n", n, ".p", p, ".k", k, ".*"))
  results <- list()
  iresult <- 1
  resultsfiles
  nfiles <- length(resultsfiles)
  for (ifile in 1:nfiles){
    load(resultsfiles[ifile])
    cat("file", ifile, "/", nfiles, ", containing", length(result), "runs\n")
    for (ires in 1:length(result)){
      results[[iresult]] <- result[[ires]]
      iresult <- iresult + 1
    }
  }
  return(results)
}
#

designs <- c(0,1)
SNRs <- c(0.5, 1, 2)
ns <- c(500,1000)
ps <- c(1000, 5000)
k <- 0
df_ <- data.frame()
for (design in designs){
  for (SNR in SNRs){
    cat("SNR=", SNR, "\n")
    for (n in ns){
      for (p in ps){
        cat("p=", p, "\n")
        results <- get_results(design, SNR, n, p, k)
        meetings <- sapply(results, function(x) x$meetingtime)
        df_ <- rbind(df_, data.frame(design = design, SNR = SNR, n = n, p = p, meanmeeting = mean(meetings), sdmeeting = sd(meetings), nrep = length(meetings)))
      }
    }
  }
}

df_$sderror <- df_$sdmeeting/sqrt(df_$nrep)

table.independent <- df_ %>% filter(design == 0) %>% select(SNR,n,p,meanmeeting) %>% arrange(n,SNR,p) %>%
  spread(SNR, meanmeeting) %>% setNames(c("n", "p", "SNR = 0.5", "SNR = 1", "SNR = 2"))

library(xtable)
cap <- "bla"
formatted.df <- xtable(table.independent, digits = 0, caption = cap)
formatted.df


table.correlated <- df_ %>% filter(design == 1) %>% select(SNR,n,p,meanmeeting) %>% arrange(n,SNR,p) %>%
  spread(SNR, meanmeeting) %>% setNames(c("n", "p", "SNR = 0.5", "SNR = 1", "SNR = 2"))
formatted.df <- xtable(table.correlated, digits = 0, caption = cap)
formatted.df

