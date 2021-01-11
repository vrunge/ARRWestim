# ARRWestim : An example


Installation of the package from github
```r
devtools::install_github("vrunge/ARRWestim")
```

## Data generation

We chose parameters
```r
n <- 1000
phi <- 0.2
sdEta <- 0.5
sdNu <- 0.2
```

And generate a time series of length `n` using the `dataARRW` function
```r
data <- dataRWAR(N = n, sdEta = sdEta, sdNu = sdNu, type = "rand1",  nbSeg = 10, seed = 8)
```

We can plot the time series
```r
plotRWAR(data)
```

and the diff-1 time-series and show the changepoints
```r
plotRWARdiff(data)
```

## Variance analysis

We robustly estimate the variances of different lag-k data
```r
nb <- 10
v <- estimVar(data$y, nbK = nb)
```

and find all parameters (the AR and RW variances and phi AR(1) autocorrelation parameter)
```r
bestParameters(data$y, nbK = nb)
```

##  Least square comparison

This last function plots the vectors to compare in the least-square criterion

```r
plotVarVarEstim(v, sdEta, sdNu, phi, nbK = nb)
```
