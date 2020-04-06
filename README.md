# ARRWestim : An example

### To install the package from github
```r
devtools::install_github("vrunge/ARRWestim")
```

### chose parameters
```r
n <- 10000
phi <- 0.2
sdEta2 <- 10
sdNu2 <- 50
```

### GENERATE DATA

```r
y <- dataARRW(n = n, sdEta2 = sdEta2, sdNu2 = sdNu2, phi = phi, poisParam = 0.001, meanGap = 2)
```

###  plot the time-series
```r
plotARRW(y)
```
###  plot the diff-1 time-series and show changepoints
```r
plotARRWdiff(y)
```

## ESTIM VARIANCES
```r
nb <- 10
v <- estimVar(y$y, nbK = nb)
```

#### find all parameters (the AR and RW variances and phi AR(1) autocorrelation parameter)
```r
bestParameters(v = v, nbK = nb)
```

###  compare data for the least-square criterion
```r
plotVarVarEstim(v, sdEta2, sdNu2, phi, nbK = nb)
```
