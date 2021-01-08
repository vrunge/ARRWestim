#' Variance estimation for diff k operators
#' @description Estimation of the variances for the diff k operator k = 1 to nbK
#' @param y A time-series obtained by the dataRWAR function
#' @param nbK number of diff k elements to consider
#' @return the vector varEst of estimated variances
#' @examples
#' estimVar(dataRWAR(1000, sdEta2 = 0.1, sdNu2 = 0.1, phi = 0.3, type = "rand1",  nbSeg = 10)$y)
estimVar <- function(y, nbK = 10)
{
  n <- length(y)
  varEst <- rep(0, nbK)
  for(k in 1:nbK)
  {
    z <- y[(k+1):n] - y[1:(n-k)]
    varEst[k] <- mad(z)^2 #k*sdEta2 + 2*((1-phi^k)/(1-phi^2))*sdNu2
  }
  return(varEst)
}


#' L2 error estimation
#' @description the least-square value
#' @param v the estimated variances of the diff k operator
#' @param sdEta2 the Random Walk variance
#' @param sdNu2 the AR(1) variance
#' @param phi the autocorrelative AR(1) parameter
#' @return the value of the sum of squares
cost <- function(v, sdEta2, sdNu2, phi)
{
  return(sum(sapply(1:length(v), function(k){(k*sdEta2 + 2*((1-phi^k)/(1-phi^2))*sdNu2 - v[k])^2})))
}


#' RW and AR(1) variance estimations with fixed AR(1) parameter
#' @description Evaluation of the variances Eta2 and Nu2
#' @param v the estimated variances of the diff k operator
#' @param phi the autocorrelative AR(1) parameter
#' @return a list with an estimation of the variances Eta2 and Nu2
evalEtaNu <- function(v, phi)
{
  nbK <- length(v)
  a1 <- sum((1:nbK)^2)
  b1 <- 2*sum((1:nbK)*(1-phi^(1:nbK))/(1-phi^2))
  a2 <- b1/2
  b2 <- 2*sum(((1-phi^(1:nbK))/(1-phi^2))^2)
  c1 <- sum((1:nbK)*v)
  c2 <- sum(((1-phi^(1:nbK))/(1-phi^2))*v)
  det <- a1*b2 - a2*b1
  myEta2 <- (b2*c1 - b1*c2)/det
  myNu2 <- (-a2*c1 + a1*c2)/det

  if(myEta2 < 0 || myNu2 <0) #the KKT condition with constraints phi with myEta2 >= 0 and myNu2 >=0
  {
    if(myEta2 < 0)
    {
      myEta2 = 0
      myNu2 = c2/b2
    }
    if(myNu2 < 0)
    {
      myEta2 = c1/a1
      myNu2 = 0
    }
  }
  return(list(Eta2 = myEta2, Nu2 = myNu2))
}


#' bestParameters
#' @description iteration of the least square criterion for a grid of the phi parameter
#' @param y A time-series obtained by the dataRWAR function
#' @param nbK number of diff k elements to consider
#' @return a list with an estimation of the best parameters for Eta2, Nu2 and phi
#' bestParameters(dataRWAR(10000, sdEta2 = 0.2, sdNu2 = 0.1, phi = 0.3, type = "rand1",  nbSeg = 10,seed = sample(10000,1))$y)
bestParameters <- function(y, nbK = 10)
{
  costall <- rep(0,100)
  v <- estimVar(y, nbK = nbK) #using function estimVar
  for(i in 1:100)
  {
    e <- evalEtaNu(v, (i-1)/100)  #using function evalEtaNu
    costall[i] <- cost(v, e$Eta2, e$Nu2, (i-1)/100)  #using cost function
  }
  argmin <- which.min(costall)
  e <- evalEtaNu(v, (argmin-1)/100)
  return(list(costall = costall,Eta2Opt = e$Eta2, Nu2Opt = e$Nu2, argmin = (argmin-1)/100))
}
