

#' one.simu
#' @description Perform the analysis form data generation to parameter estimation
#' @param i index of the simulation
#' @param N number of data points
#' @param sdEta2 variance in Random Walk
#' @param sdNu2 variance in AR(1)
#' @param phi AR(1) autocorrelation parameter
#' @param type possible scenarios for the jump structure
#' @param nbSeg number of segment
#' @param jumpSize max jump size
#' @param nbK number of diff k elements to consider
#' @return a dataframe with the inital parameters and the estimated ones
one.simu <- function(i, N = 10^5, sdEta2 = 0.4, sdNu2 = 0.3, phi = 0.2,  type = "rand1", nbSeg = 10, jumpSize = 2, nbK = 10)
{
  y <- dataRWAR(N = N,
                sdEta2 = sdEta2, sdNu2 = sdNu2, phi = phi,
                type = type,
                nbSeg = nbSeg, jumpSize = jumpSize,
                seed = sample(1e6,1))
  res <- bestParameters(y$y, nbK = nbK)

  #ind, n, phi, sdEta2, sdNu2, nbK, poiP, meanGap, method, mse, beta, K
  df <- data.frame(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0) ,stringsAsFactors = FALSE)
  colnames(df) <- c("index", "n", "sdEta2", "sdNu2", "phi", "nbK",
                    "nbSeg", "jumpSize", "sdEta2Est%", "sdNu2Est%", "phiEst_error")
  df[1,] <- c(i,N, sdEta2, sdNu2, phi, nbK, nbSeg, jumpSize,
              (res$Eta2Opt - sdEta2)/sdEta2, (res$Nu2Opt-sdNu2)/sdNu2, res$argmin - phi)
  return(df)
}
