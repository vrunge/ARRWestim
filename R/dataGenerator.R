#' Signal scenario
#'
#' @description Generate a signal scenario
#' @param N number of data points
#' @param type possible scenarios for the jump structure
#' @param nbSeg number of segments
#' @param jumpSize max size of the jumps
#' @param seed random number generator seed (default is 42)
scenarioGenerator <- function(N, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 1, seed = 42)
{
  #segment length
  set.seed(seed)
  type <- match.arg(type)
  if (type == "rand1")
  {
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
    else if(s < N)
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
  }

  switch(
    type,
    none = rep(0, N),
    up = unlist(lapply(0:(nbSeg-1), function (k) rep(k * jumpSize, N * 1 / nbSeg))),
    updown = unlist(lapply(0:(nbSeg-1), function(k) rep((k %% 2) * jumpSize, N * 1 / nbSeg))),
    rand1 = unlist(sapply(1:nbSeg, function(i) rep(rand1Jump[i] * jumpSize, r1[i])))
  )
}


#' Generating data function with the deCAFS model
#'
#'
#' @description Generate a Realization from the deCAFS model (see references for further details).
#' \deqn{y_t = \mu_t + \epsilon_t}
#' where
#' \deqn{\mu_t = \mu_{t-1} + \eta_t + \delta_t, \quad \eta_t \sim N(0, \sigma_\eta^2), \ \delta_t \ \in R}
#' and
#' \deqn{\epsilon_t = \phi \epsilon_{t-1} + \nu_t \quad \nu_t \sim N(0, \sigma_\nu^2)}
#' @param N number of data points
#' @param sdEta standard deviation in Random Walk
#' @param sdNu standard deviation in AR(1)
#' @param phi AR(1) autocorrelation parameter
#' @param type possible scenarios for the jump structure
#' @param nbSeg number of segments
#' @param jumpSize max size of the jumps
#' @param seed random number generator seed (default is 42)
#'
#' @return A list containing:
#' \describe{
#' \item{\code{y}}{the data sequence,}
#' \item{\code{signal}}{the underlying signal without the superimposed AR(1) noise,}
#' \item{\code{changepoints}}{the changepoint indices}
#' }
#' @examples
#' myData <- dataRWAR(1000, sdEta = 0.1, sdNu = 0.1, type = "rand1",  nbSeg = 10, seed = 86)

dataRWAR <- function(N = 1e3,
                     sdEta = 1, sdNu = 1, phi = 0.5,
                     type = c("none", "up", "updown", "rand1"),
                     nbSeg = 20, jumpSize = 1, seed = 42)
{
  model <- scenarioGenerator(N, type = type, nbSeg = nbSeg, jumpSize = jumpSize, seed = seed)
  #Randow walk element
  RW <- cumsum(rnorm(n = N, 0, sdEta))
  y <- model + RW

  #AR(1) noise
  e <- rnorm(1, 0, sdNu/sqrt(1 - phi^2))
  for(i in 1:N)
  {
    y[i] <- y[i] + e
    e <- phi*e + rnorm(n = 1, mean = 0, sd = sdNu)
  }
  return(list(y = y, signal = model + RW, changepoints = which(diff(model) != 0)))
}

