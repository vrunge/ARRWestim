#' dataARRW
#' @description Generate data with the ARRW model
#' @param n number of data points
#' @param sdEta2 variance in Random Walk
#' @param sdNu2 variance in AR(1)
#' @param phi AR(1) autocorrelation parameter
#' @param poisParam a probability in Poisson distribution for changepoint density
#' @param meanGap mean gap of the jumps
#' @return a vector of simulated data with y = the data, mu = the signal, changepoints = changepoint indices
#' @examples
#' myData <- dataARRW(100)
dataARRW <- function(n = 1e3, sdEta2 = 1, sdNu2 = 1, phi = 0.9, poisParam = 0, meanGap = 2)
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
