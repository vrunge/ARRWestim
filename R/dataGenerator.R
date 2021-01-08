#' dataARRW
#' @description Generate data with the ARRW model
#' @param model the initial model for change-point position and jump size
#' @param sdEta2 variance in Random Walk
#' @param sdNu2 variance in AR(1)
#' @param phi AR(1) autocorrelation parameter
#' @return a vector of simulated data with y = the data, mu = the signal, changepoints = changepoint indices
#' @examples
#' myData <- dataARRW(scenarioGenerator(100, "rand1", nbSeg = 20))
dataARRW <- function(model, sdEta2 = 1, sdNu2 = 1, phi = 0.9)
{
  n <- length(model)
  RW <- cumsum(rnorm(n, 0, sqrt(sdEta2)))

  y <- model + RW
  e <- rnorm(1, 0, sqrt(sdNu2)/sqrt(1 - phi^2))
  for(i in 1:n)
  {
    y[i] <- y[i] + e
    e <- phi*e + rnorm(n = 1, mean = 0, sd = sqrt(sdNu2))
  }
  return(list(y = y, mu = model + RW, changepoints = which(diff(model) != 0)))
}



#' scenarioGenerator
#' @description Generate model for ARRW
#' @param N number of data points
#' @param type 4 possible scenario
#' @param nbSeg number of segments
#' @param jumpSize Max size of the jumps
#' @param seed default seed is 42
#' @return a vector of piecewise constant segment
#' @examples
#' s <- scenarioGenerator(100, "rand1", nbSeg = 20)
scenarioGenerator <- function(N, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 1, seed = 42)
{
  #segment length
  set.seed(seed)
  rand1CP <- rpois(nbSeg, lambda = 10)
  r1 <- pmax(round(rand1CP * N / sum(rand1CP)), 1) #normalisation (delete 0 values)
  s <- sum(r1)
  if(s > N)
  {
    while(sum(r1) > N)
    {
      p <- sample(x = nbSeg, size = 1)
      if(r1[p]> 1){r1[p] <- r1[p] - 1}
    }
  }

  if(s < N)
  {
    for(i in 1:(N-s))
    {
      p <- sample(x = nbSeg, size = 1)
      r1[p] <- r1[p] + 1
    }
  }

  #jump intensity
  set.seed(seed + 1)
  rand1Jump <- runif(nbSeg, min = -1, max = 1)

  type <- match.arg(type)
  switch(
    type,
    none = rep(0, N),
    up = unlist(lapply(0:(nbSeg-1), function (k) rep(k * jumpSize, N * 1 / nbSeg))),
    updown = unlist(lapply(0:(nbSeg-1), function(k) rep((k %% 2) * jumpSize, N * 1 / nbSeg))),
    rand1 = unlist(sapply(1:nbSeg, function(i) rep(rand1Jump[i] * jumpSize, r1[i])))
  )
}

