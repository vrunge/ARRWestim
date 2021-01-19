rm(list=ls())

# LIBRARIES
library(ARRWestim)
library(parallel)
library(fields)

###
### NB cores for parallel computing
###
cores <- detectCores()
cores <- 45
#cores <- 8 ###CHANGE

###
### Simulations parameters
###
#nbSimu <- 1200
nbSimu <- 1000 ##CHANGE
nbPhi <- 19 #step size = 0.05 in phi
nbOmega2 <- 40
nbK <- 10


###
### phi and omega2 = grid for simulations (omega2 = (sd_eta/sd_nu)^2)
### sd_nu fixed to 1
###

phi <- seq(from = 0, to = 0.9, length.out = nbPhi)

omega2 <- exp(seq(from = -log(3), to = log(8), length.out = nbOmega2))
diffO2 <- diff(log(omega2))[1]
logOmega2 <- c(log(omega2)[1]-diffO2, log(omega2))
omega2 <- c(0, omega2)
nbOmega2 <- nbOmega2 +1

#omega2 regular on log scale


###
### label and positions for y-axis image plot
###
eps <- -2*10^(-15)
nb0.5 <- which(omega2-0.5 >= eps)[1]
nb1 <- which(omega2-1 >= eps)[1]
nb4 <- which(omega2-4 >= eps)[1]
nb8 <- which(omega2-8 >= eps)[1]
myscale <- c(0,0.5,1,4,8)
positions <- c(logOmega2[1], logOmega2[nb0.5], logOmega2[nb1], logOmega2[nb4], logOmega2[nb8])


######################################################################
# color scale function

colScale <- function(min, max, nb, epsilon)
{
  if(min*max < 0)
  {
    colorTable <- designer.colors(3*nb-1, c( "blue","white", "red"))
    M <- max(abs(min), abs(max))
    brks<- c(seq(-M, -epsilon, length.out = nb+1)[-(nb+1)],
             seq(-epsilon, epsilon,length.out = nb+1)[-(nb+1)],
             seq( epsilon, M,length.out = nb))
    if(-M < min)
    {
      rk <- sum(brks < min)
      colorTable <- colorTable[rk:(3*nb-1)]
      brks <- brks[rk:(3*nb)]
    }
    else
    {
      rk <- sum(brks > max)
      colorTable <- colorTable[1:(3*nb-1-rk)]
      brks <- brks[1:(3*nb-rk)]
    }
  }
  if(min*max >= 0)
  {
    if(min < 0)
    {
      colorTable <- designer.colors(nb, c( "blue","white"))
      brks<- seq(min, 0, length.out = nb)
    }
    else
    {
      colorTable <- designer.colors(nb, c( "white","red"))
      brks<- seq(0, max, length.out = nb)
    }

  }
  return(list(col = colorTable, breaks = brks))
}


grayScale <- function(max, nb)
{
  colorTable <- designer.colors(nb-1, c("white", "black"))
  brks<- seq(0, max, length.out = nb)
  return(list(col = colorTable, breaks = brks))
}


########### ########### ########### ###########

res1 <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res1 <- c(res1, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "none",
                             nbSeg = 1,
                             jumpSize = 0,
                             nbK = nbK,
                             varType = "S",
                             mc.cores = cores)) ## mc.cores = 8
  }
}





res2 <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res2 <- c(res2, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 5,
                             jumpSize = 10,
                             nbK = nbK,
                             varType = "S",
                             mc.cores = cores)) ## mc.cores = 8
  }
}


res3 <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res3 <- c(res3, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 25,
                             jumpSize = 10,
                             nbK = nbK,
                             varType = "S",
                             mc.cores = cores)) ## mc.cores = 8
  }
}

########### ########### ########### ###########

res1MAD <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res1MAD <- c(res1MAD, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "none",
                             nbSeg = 1,
                             jumpSize = 0,
                             nbK = nbK,
                             varType = "MAD",
                             mc.cores = cores)) ## mc.cores = 8
  }
}





res2MAD <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res2MAD <- c(res2MAD, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 5,
                             jumpSize = 10,
                             nbK = nbK,
                             varType = "MAD",
                             mc.cores = cores)) ## mc.cores = 8
  }
}


res3MAD <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res3MAD <- c(res3MAD, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 25,
                             jumpSize = 10,
                             nbK = nbK,
                             varType = "MAD",
                             mc.cores = cores)) ## mc.cores = 8
  }
}



df1_S <- do.call(rbind, res1)
df2_S <- do.call(rbind, res2)
df3_S <- do.call(rbind, res3)

df1_MAD <- do.call(rbind, res1MAD)
df2_MAD <- do.call(rbind, res2MAD)
df3_MAD <- do.call(rbind, res3MAD)
