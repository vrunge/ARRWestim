# ARRWestim : An example


Installation of the package from github
```r
devtools::install_github("vrunge/ARRWestim")
```

## Data generation

We chose parameters
```r
n <- 10000
phi <- 0.2
sdEta2 <- 10
sdNu2 <- 50
```

And generate a time series of length `n` using the `dataARRW` function
```r
y <- dataARRW(n = n, sdEta2 = sdEta2, sdNu2 = sdNu2, phi = phi, poisParam = 0.001, meanGap = 2)
```

We can plot the time series
```r
plotARRW(y)
```

and the diff-1 time-series and show the changepoints
```r
plotARRWdiff(y)
```

## Variance analsis

We robustly estimate the variances of different lag-k data
```r
nb <- 10
v <- estimVar(y$y, nbK = nb)
```

and find all parameters (the AR and RW variances and phi AR(1) autocorrelation parameter)
```r
bestParameters(v = v, nbK = nb)
```

##  Least square comparison

This last function plots the vectors to compare in the least-square criterion

```r
plotVarVarEstim(v, sdEta2, sdNu2, phi, nbK = nb)
```
