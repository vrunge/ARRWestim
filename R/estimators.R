
#' estimVar
#' @description Estimation of the variances for the diff k operator k = 1 to nbK
#' @param y An object created by the dataARRW function
#' @param nbK number of diff k elements to consider
#' @return a list with a varEst element (= a vector)
estimVar <- function(y, nbK = 10)
{
  n <- length(y)
  varEst <- rep(0, nbK)
  for(k in 1:nbK)
  {
    z <- y[(k+1):n] - y[1:(n-k)]
    varEst[k] <- mad(z)^2 #k*sdEta2 + 2*((1-phi^k)/(1-phi^2))*sdNu2
  }
  return(list(varEst = varEst))
}



#' cost
#' @description the least-square value
#' @param v the estimated variances of the diff k operator
#' @param sdEta2 the Random Walk variance
#' @param sdNu2 the AR(1) variance
#' @param phi the autocorrelative AR(1) parameter
#' @param nbK number of diff k elements to consider
#' @return the value of the sum of squares
cost <- function(v, sdEta2, sdNu2, phi, nbK = 10)
{
  v <- unname(unlist(v))
  return(sum(sapply(1:nbK, function(k){(k*sdEta2 + 2*((1-phi^k)/(1-phi^2))*sdNu2 - v[k])^2})))
}


#' evalEtaNu
#' @description Evaluation of the variances Eta2 and Nu2
#' @param v the estimated variances of the diff k operator
#' @param phi the autocorrelative AR(1) parameter
#' @param nbK number of diff k elements to consider
#' @return a list with an estimation of the variances Eta2 and Nu2
evalEtaNu <- function(v, phi, nbK = 10)
{
  v <- unname(unlist(v))
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
#' @param v the estimated variances of the diff k operator
#' @param nbK number of diff k elements to consider
#' @return a list with an estimation of the best parameters for Eta2, Nu2 and phi
bestParameters <- function(v, nbK = 10)
{
  costall <- rep(0,100)
  for(i in 1:100)
  {
    e <- evalEtaNu(v, (i-1)/100, nbK = nbK)
    costall[i] <- cost(v, e$Eta2, e$Nu2, (i-1)/100, nbK)
  }
  argmin <- which.min(costall)
  e <- evalEtaNu(v, (argmin-1)/100, nbK = nbK)
  return(list(Eta2Opt = e$Eta2, Nu2Opt = e$Nu2, argmin = (argmin-1)/100))
}



