#' l2Threshold
#' @description Simple linear threshold changepoint detection method
#' @param data a time-series
#' @param beta the penalty value
#' @param lambda 1/variance2
#' @return a list of changepoints with their jump sizes
l2Threshold <- function(data, beta, lambda)
{
  n <- length(data)
  test <- lambda * (data[-1] - data[-n])^2
  z <- (test > beta) * (data[-1] - data[-n])
  chgpt <- which(z != 0)
  jump <- z[chgpt]
  return(list(changepoints = c(chgpt+1, n), jumpSize = jump))
}
