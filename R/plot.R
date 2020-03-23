###
###  plotARRW
###

plotARRW <- function(y)
{
  ylim <- c(min(y$mu,y$y), max(y$mu,y$y))
  plot(y$y, col = 1, ylim = ylim, pch = '+', ylab = "")# y in blue
  par(new= TRUE)
  plot(y$mu, col = 2, ylim = ylim, pch = '+', ylab = "") # true signal in red
  abline(v = y$changepoints)
}


###
###  plotARRWdiff
###

plotARRWdiff <- function(y)
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




plotVarVarEstim <- function(v, sdEta2, sdNu2, phi, nbK = 10)
{
  v <- v$varEst
  #### ESTIM var
  vari <- rep(0,nbK)
  for(k in 1:nbK)
  {
    vari[k] <- k*sdEta2 + 2*((1-phi^k)/(1-phi^2))*sdNu2
  }

  ylim <- c(min(vari,v), max(vari,v))
  plot(vari, ylim = c(ylim), col = 1)
  par(new = TRUE)
  plot(v, ylim = c(ylim), col = 2) # red = estimation
}
