

estimVar <- function(y, nbK = 10)
{
  varEst <- rep(0, nbK)
  for(k in 1:nbK)
  {
    z <- y[(k+1):n] - y[1:(n-k)]
    varEst[k] <- mad(z)^2 #k*sdEta2 + 2*((1-phi^k)/(1-phi^2))*sdNu2
  }
  return(list(varEst = varEst))
}



###
###  least square cost function
###

cost <- function(v, sdEta2, sdNu2, phi, nbK = 10)
{
  v <- unname(unlist(v))
  return(sum(sapply(1:nbK, function(k){(k*sdEta2 + 2*((1-phi^k)/(1-phi^2))*sdNu2 - v[k])^2})))
}


###
###  evalEtaNu with a given phi
###

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
  return(list(Eta2 = myEta2, Nu2 = myNu2))
}




###
###  bestParameters
###

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



