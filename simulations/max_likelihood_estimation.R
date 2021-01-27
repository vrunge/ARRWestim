rm(list=ls())

# LIBRARIES
#devtools::install_github("gtromano/DeCAFS")
#devtools::install_github("vrunge/ARRWestim")
library(DeCAFS)
library(ARRWestim)


N <- 5e3
sdEta <- 0.1
sdNu <- 0.3
phi <- 0.3
type <- "rand1"
nbSeg <- 10
jumpSize <- 2
nbK <- 10
varType <- "MAD"

### DATA GENERATION

y <- ARRWestim::dataRWAR(N = N,
                         sdEta = sdEta, sdNu = sdNu, phi = phi,
                         type = type,
                         nbSeg = nbSeg, jumpSize = jumpSize,
                         seed = sample(1e6,1))

### PARAMETER ESTIMATION (our method in the paper)
res <- bestParameters(y$y, nbK = nbK, type = varType)
res

### DeCAFS
deca <- DeCAFS::DeCAFS(y$y, 2*log(N),
                       list(sdEta = res$EtaOpt, sdNu = res$NuOpt, phi = res$argmin))

### MLE for PHI with the estimated signal
Z <- y$y - deca$signal
MLE_phi <- sum(Z[-1]*Z[-N])/(sum(Z[c(-1,-N)]*Z[c(-1,-N)]))

#comparison
MLE_phi
res$argmin

### MLE phi with the true signal y$signal
Z <- y$y - y$signal
MLE_phi_trueSignal <- sum(Z[-1]*Z[-N])/(sum(Z[c(-1,-N)]*Z[c(-1,-N)]))
MLE_phi_trueSignal
