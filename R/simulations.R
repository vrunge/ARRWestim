

#' one.simu
#' @description Perform the analysis form data generation to parameter estimation
#' @param i index of the simulation
#' @param n number of data points
#' @param sdEta2 variance in Random Walk
#' @param sdNu2 variance in AR(1)
#' @param phi AR(1) autocorrelation parameter
#' @param poisParam a probability in Poisson distribution for changepoint density
#' @param meanGap mean gap of the jumps
#' @param nbK number of diff k elements to consider
#' @return a dataframe with the inital parameters and the estimated ones
one.simu <- function(i, n = 10^5, sdEta2 = 1, sdNu2 = 1, phi = 0.5, poisParam = 0, meanGap = 2, nbK)
{
  y <- dataARRW(n = n, sdEta2 = sdEta2, sdNu2 = sdNu2,  phi = phi, poisParam, meanGap)
  v <- estimVar(y = y$y, nbK = nbK) ##vector of estimated variances
  res <- bestParameters(v = v, nbK = nbK)

  #ind, n, phi, sdEta2, sdNu2, nbK, poiP, meanGap, method, mse, beta, K
  df <- data.frame(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0) ,stringsAsFactors = FALSE)
  colnames(df) <- c("index", "n", "sdEta2", "sdNu2", "phi", "nbK", "poisParam", "meanGap", "sdEta2Est%", "sdNu2Est%", "phiEst")
  df[1,] <- c(i,n, sdEta2, sdNu2, phi, nbK, poisParam, meanGap,
              (res$Eta2Opt - sdEta2)/sdEta2, (res$Nu2Opt-sdNu2)/sdNu2, res$argmin - phi)
  return(df)
}
