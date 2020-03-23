# ARRWestim : an example

### chose parameters
n <- 10000
phi <- 0.2
sdEta2 <- 10
sdNu2 <- 50

### GENERATE DATA
y <- dataARRW(n = n, poisParam = 0.001, meanGap = 2, phi = phi, sdEta2 = sdEta2, sdNu2 = sdNu2)

###  plot the time-series
plotARRW(y)

###  plot the diff-1 time-series and show changepoints
plotARRWdiff(y)

## ESTIM VARIANCES
nb <- 10
v <- estimVar(y$y, nbK = nb)

#### find all parameters
bestParameters(v = v, nbK = nb)


###  compare data for the least-square criterion
plotVarVarEstim(v, sdEta2, sdNu2, phi, nbK = nbK)

