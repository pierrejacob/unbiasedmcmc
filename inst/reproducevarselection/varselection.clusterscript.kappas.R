# time allocated to this
TIME <- 2*60*60

library(unbiasedmcmc)

library(parallel)
# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# igrid <- 1
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

n <- 500
p <- 1000
SNR <- 1
#
# load data
load(paste0("varselection.dataSNR", SNR, ".RData"))
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

ikappa <- 1 + (igrid-1) %% nkappas
kappa <- kappas[ikappa]

k <- 75000
m <- 150000

print(ikappa)

# load model
vs <- get_variableselection(Y,X,g,kappa,s0,proportion_singleflip)
prior <- vs$prior
marginal_likelihood <- vs$marginal_likelihood
rinit <- vs$rinit
single_kernel <- vs$single_kernel
coupled_kernel <- vs$coupled_kernel

unbiasedestimator <- function(single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 1, max_iterations = Inf, totalduration = Inf){
  ptm <- proc.time()
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  current_pdf1 <- marginal_likelihood(chain_state1) + prior(chain_state1)
  current_pdf2 <- marginal_likelihood(chain_state2) + prior(chain_state2)
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(chain_state1)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  sres1 <- single_kernel(chain_state1, current_pdf1)
  chain_state1 <- sres1$state
  current_pdf1 <- sres1$pdf
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(chain_state1)
  }
  # Check if time is up already
  elapsedtime <- as.numeric((proc.time() - ptm)[3])
  if (elapsedtime > totalduration){
    # time's up
    return(list(finished = FALSE, message = "interrupted because time's up"))
  }
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    elapsedtime <- as.numeric((proc.time() - ptm)[3])
    if (elapsedtime > totalduration){
      # time's up
      return(list(finished = FALSE, message = "interrupted because time's up"))
    }
    iter <- iter + 1
    if (meet){
      sres1 <- single_kernel(chain_state1, current_pdf1)
      chain_state1 <- sres1$state
      current_pdf1 <- sres1$pdf
      chain_state2 <- chain_state1
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(chain_state1)
      }
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
      chain_state1 <- res_coupled_kernel$state1
      current_pdf1 <- res_coupled_kernel$pdf1
      chain_state2 <- res_coupled_kernel$state2
      current_pdf2 <- res_coupled_kernel$pdf2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
      }
    }
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction
  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}


ptm_irep <- proc.time()
last_start <- ptm_irep
nsamples_ <- 0
durations_ <- c()
elapsedtime <- 0
timeleft <- TIME
results_ <- list()
while(elapsedtime < TIME){
  if (nsamples_ == 0){ # then produce a sample
    res <- try(unbiasedestimator(single_kernel, coupled_kernel, rinit, h = function(x) x, k = k, m = m))
    if (inherits(res, "try-error")){
      res <- list(finished = FALSE, meeting = Inf, uestimator = rep(0, p))
    }
  } else { # then produce a sample if time permits
    res <- try(unbiasedestimator(single_kernel, coupled_kernel, rinit, h = function(x) x, k = k, m = m, totalduration = timeleft))
    if (inherits(res, "try-error")){
      res <- list(finished = FALSE, meeting = Inf, uestimator = rep(0, p))
    }
  }
  elapsedtime <- as.numeric((proc.time() - ptm_irep)[3])
  timeleft <- TIME - elapsedtime
  duration_ <- as.numeric((proc.time() - last_start)[3])
  last_start <- proc.time()
  if (res$finished){
    nsamples_ <- nsamples_ + 1
    results_[[nsamples_]] <- res
    durations_ <- c(durations_, duration_)
  }
  save(results_, nsamples_, durations_, file = paste0("output/varselection.ikappa", ikappa, ".ID", igrid, ".RData"))
}

