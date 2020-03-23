
l2Threshold <- function(data, beta, lambda)
{
  n <- length(data)
  test <- lambda * (data[-1] - data[-n])^2
  z <- (test > beta) * (data[-1] - data[-n])

  chgpt <- which(z != 0)
  jump <- z[chgpt]

  return(list(changepoints = c(chgpt+1, n), jumpSize = jump))
}
