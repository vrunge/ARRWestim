###
###  dataARRW
###

dataARRW <- function(n = 1e3, poisParam = 0, meanGap = 10, phi = .9, sdEta2 = 1, sdNu2 = 1)
{
  changepoints <- rpois(n, poisParam)
  changes <- cumsum(sample(c(-1, 1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap))
  RW <- cumsum(rnorm(n, 0, sqrt(sdEta2)))
  res <- changes + RW

  e <- rnorm(1, 0, sqrt(sdNu2)/sqrt(1 - phi^2))
  for(i in 1:n)
  {
    res[i] <- res[i] + e
    e <- phi*e + rnorm(n = 1, mean = 0, sd = sqrt(sdNu2))
  }

  return(list(y = res, mu = changes + RW, changepoints = which(changepoints > 0)))
}
