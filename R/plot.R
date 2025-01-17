#' plotARRW
#' @description plot signal, y datapoints and changepoints on a same graph
#' @param y An object created by the dataARRW function
plotRWAR <- function(y)
{
  ylim <- c(min(y$signal,y$y), max(y$signal,y$y))
  plot(y$y, col = 1, ylim = ylim, pch = '+', ylab = "")# y in blue
  par(new= TRUE)
  plot(y$signal, col = 2, ylim = ylim, pch = '+', ylab = "") # true signal in red
  abline(v = y$changepoints)
}

#' plotRWARdiff
#' @description plot y_{t+1} - y_t signal and changepoint locations to identify obvious changepoints
#' @param y An object created by the dataARRW function
plotRWARdiff <- function(y)
{
  z <- diff(c(0,y$y))
  plot(z,xlim = c(1,length(z)))
  if(length(y$changepoints) > 0)
  {
    par(new = TRUE)
    plot(y$changepoints,z[y$changepoints], xlim = c(1,length(z)), col = 2)
    abline(v = y$changepoints)
  }
}



#' plotVarVarEstim
#' @description plot the estimated variances v_k against the true variances for the diff k operator (y_{t+k} - y_t) for k = 1 to nbK
#' @param v the estimated variances of the diff k operator
#' @param sdEta standard deviation in Random Walk
#' @param sdNu standard deviation in AR(1)
#' @param phi the autocorrelative AR(1) parameter
#' @param nbK number of diff k elements to consider
plotVarVarEstim <- function(v, sdEta, sdNu, phi, nbK = 10)
{
  #### ESTIM var
  vari <- rep(0,nbK)
  for(k in 1:nbK)
  {
    vari[k] <- k*sdEta^2 + 2*((1-phi^k)/(1-phi^2))*sdNu^2
  }

  ylim <- c(min(vari,v), max(vari,v))
  plot(vari, ylim = c(ylim), col = 1)
  par(new = TRUE)
  plot(v, ylim = c(ylim), col = 2) # red = estimation
}
