# send to Jay Emerson? Louis Aslett?

dimstate = 100
nmcmc = 1e4
chain = matrix(0, nrow = nmcmc, ncol = dimstate)
for (imcmc in 1:nmcmc){
  if (imcmc == nrow(chain)){
  }
  x = rnorm(dimstate, mean = 0, sd = 1)
  chain[imcmc,] = x
}

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

dimstate = 100
nmcmc = 1e4
chain = matrix(0, nrow = nmcmc, ncol = dimstate)
for (imcmc in 1:nmcmc){
  if (imcmc == dim(chain)[1]){
  }
  x = rnorm(dimstate, mean = 0, sd = 1)
  chain[imcmc,] = x
}


x <- matrix(0, nrow=1e5, ncol=100) # matrix has ref count 1
x[1,1] <- 1 # ref count is 1, so write with no copy
nrow(x) # ref count is 2 even though nothing was touched
x[1,1] <- 1 # ref count still 2, so R copies before writing first element. Now the ref count drops to 1 again
x[2,2] <- 1 # this writes without a copy as ref count got reset on last line
nrow(x) # ref count jumps
x[3,3] <- 1 # copy invoked again! Aaaargh!
