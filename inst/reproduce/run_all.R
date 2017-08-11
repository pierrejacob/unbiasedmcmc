# this folder
scriptfolder <- "~/Dropbox/PolyaGamma/code/debiasedmcmc/inst/reproduce/"
figurefolder <- "~/Dropbox/PolyaGammaResults/reproduce/"
resultsfolder <- "~/Dropbox/PolyaGammaResults/reproduce/"

# note: we could automatize by running all R scripts finishing in run.R
# note: we should remove "setwd" from the scripts themselves, and handle it from here

setwd(scriptfolder)
source("mixture.run.R")
setwd(scriptfolder)
source("plummer.cut.run.R")
setwd(scriptfolder)
source("pumpfailures.run.R")
setwd(scriptfolder)
source("baseball.run.R")
setwd(scriptfolder)
# ...
setwd(scriptfolder)
# ...
setwd(scriptfolder)
# ...
setwd(scriptfolder)
