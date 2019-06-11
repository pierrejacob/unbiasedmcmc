# This script illustrates a weird phenomenon when using
# nrow/ncol to get a matrix' dimension instead of dim.

# The mystery has been elucidated thanks to Louis Aslett,
# and explanations are given in the blog post
# https://statisfaction.wordpress.com/2017/12/10/nrow-references-and-copies/

# The question is: why is the following code so slow?
dimstate = 100
nmcmc = 1e4
chain = matrix(0, nrow = nmcmc, ncol = dimstate)
for (imcmc in 1:nmcmc){
  if (imcmc == nrow(chain)){
  }
  x = rnorm(dimstate, mean = 0, sd = 1)
  chain[imcmc,] = x
}

# Attempts at finding the reason, identifying "nrow" as the problem
# in combination to changing the matrix chain.
dimstate = 100
nmcmc = 1e4
chain = matrix(0, nrow = nmcmc, ncol = dimstate)
for (imcmc in 1:nmcmc){
  if (imcmc == nrow(chain)){
  }
  x = rnorm(dimstate, mean = 0, sd = 1)
  # chain[imcmc,] = x
}

dimstate = 100
nmcmc = 1e4
chain = matrix(0, nrow = nmcmc, ncol = dimstate)
for (imcmc in 1:nmcmc){
  if (imcmc == nmcmc){
  }
  x = rnorm(dimstate, mean = 0, sd = 1)
  chain[imcmc,] = x
}

# illustration that dim behaves very differently compared to nrow
dimstate = 100
nmcmc = 1e4
chain = matrix(0, nrow = nmcmc, ncol = dimstate)
for (imcmc in 1:nmcmc){
  if (imcmc == dim(chain)[1]){
  }
  x = rnorm(dimstate, mean = 0, sd = 1)
  chain[imcmc,] = x
}

#
x <- matrix(0, nrow=1e5, ncol=100) # matrix has ref count 1
x[1,1] <- 1 # ref count is 1, so write with no copy
nrow(x) # ref count is 2 even though nothing was touched
x[1,1] <- 1 # ref count still 2, so R copies before writing first element. Now the ref count drops to 1 again
x[2,2] <- 1 # this writes without a copy as ref count got reset on last line
nrow(x) # ref count jumps
x[3,3] <- 1 # copy invoked again! Aaaargh!
