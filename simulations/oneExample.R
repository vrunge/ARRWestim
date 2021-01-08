### chose parameters
n <- 10000
sdEta2 <- 0.04
sdNu2 <- 0.02
phi <- 0.4

### GENERATE DATA
y <- dataRWAR(N = n, sdEta2 = sdEta2, sdNu2 = sdNu2, phi = phi, type = "rand1",  nbSeg = 20, seed = sample(1e5,1))

###  plot the time-series
plotARRW(y)

###  plot the diff-1 time-series and show changepoints
plotARRWdiff(y)

## ESTIM VARIANCES
nb <- 10
v <- estimVar(y$y, nbK = nb)
v

#### find all parameters (the AR and RW variances and phi AR(1) autocorrelation parameter)
bestParameters(y$y, nbK = nb)

###  compare data for the least-square criterion
plotVarVarEstim(v, sdEta2, sdNu2, phi, nbK = nb)

