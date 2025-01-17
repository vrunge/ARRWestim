### chose parameters
n <- 5000
sdEta <- 2
sdNu <- 1
phi <- 0.4

### GENERATE DATA
y <- dataRWAR(N = n, sdEta = sdEta, sdNu = sdNu, phi = phi,
              type = "rand1",  nbSeg = 50, jumpSize = 10,
              seed = sample(1e5,1))

###  plot the time-series
plotRWAR(y)

###  plot the diff-1 time-series and show changepoints
plotRWARdiff(y)

## ESTIM VARIANCES
nb <- 10
v1 <- estimVar(y$y, nbK = nb)
v2 <- estimVar(y$y, nbK = nb, type = "S")

#### find all parameters (the AR and RW variances and phi AR(1) autocorrelation parameter)
bestParameters(y$y, nbK = nb, type = "MAD")
bestParameters(y$y, nbK = nb, type = "S")

###  compare data for the least-square criterion
plotVarVarEstim(v1, sdEta, sdNu, phi, nbK = nb)
plotVarVarEstim(v2, sdEta, sdNu, phi, nbK = nb)
