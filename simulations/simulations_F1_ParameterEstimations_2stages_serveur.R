rm(list=ls())

# LIBRARIES
#devtools::install_github("gtromano/DeCAFS")
library(DeCAFS)
library(ARRWestim)

one.simu.2stages <- function(i, N = 5*10^3, sdEta = 0.1, sdNu = 0.3, phi = 0.2,
                     type = "rand1", nbSeg = 10,
                     jumpSize = 2, nbK = 10, varType = "MAD")
{
  y <- ARRWestim::dataRWAR(N = N,
                sdEta = sdEta, sdNu = sdNu, phi = phi,
                type = type,
                nbSeg = nbSeg, jumpSize = jumpSize,
                seed = sample(1e6,1))
  #1st estimate
  res <- bestParameters(y$y, nbK = nbK, type = varType)
  #DECAFS
  deca <- DeCAFS::DeCAFS(y$y, 2*log(N),
      list(sdEta = res$EtaOpt, sdNu = res$NuOpt, phi = res$argmin))


  chpts <- c(0, deca$changepoints, N)
  K <- length(chpts) - 1
  w <- diff(chpts)

  newRes <- matrix(0, nrow = 0, ncol = 5) #start + end + 3 parameters estimated
  for(i in 2:length(chpts))
  {
    if(w[i-1] > 50){newRes <- rbind(newRes, c(chpts[(i-1):i],0,0,0))}
  }
  if(nrow(newRes) > 0)
  {
    for(i in 1:nrow(newRes))
    {
      seg <- y$y[(newRes[i,1]+1):newRes[i,2]]
      newRes[i,3:5] <- unname(unlist(bestParameters(seg, nbK = nbK, type = varType)))
    }
  }
  finalRes <- apply(newRes, weighted.mean, newRes[,2] - newRes[,1], MARGIN = 2)
  #ind, n, phi, sdEta2, sdNu2, nbK, poiP, meanGap, method, mse, beta, K
  df <- data.frame(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0) ,stringsAsFactors = FALSE)
  colnames(df) <- c("index", "n", "sdEta", "sdNu", "phi", "nbK",
                    "nbSeg", "jumpSize", "sdEtaEst%", "sdNuEst%", "phiEst_error")
  df[1,] <- c(i,N, sdEta, sdNu, phi, nbK, nbSeg, jumpSize,
              (finalRes[3] - sdEta)/sdEta, (finalRes[4]-sdNu)/sdNu, finalRes[5] - phi)
  return(df)
}

###################

library(parallel)
library(fields)

###
### NB cores for parallel computing
###
cores <- detectCores()
cores <- 20
#cores <- 8 ###CHANGE

###
### Simulations parameters
###
#nbSimu <- 1200
nbSimu <- 1000 ##CHANGE
nbPhi <- 18 #step size = 0.05 in phi
nbOmega2 <- 40
nbK <- 10


###
### phi and omega2 = grid for simulations (omega2 = (sd_eta/sd_nu)^2)
### sd_nu fixed to 1
###

phi <- seq(from = 0, to = 0.85, length.out = nbPhi)

omega2 <- exp(seq(from = -log(12), to = log(2), length.out = nbOmega2))
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
nb4 <- which(omega2-1.5 >= eps)[1]
nb8 <- which(omega2-2 >= eps)[1]
myscale <- c(0,0.5,1,1.5,2)
positions <- c(logOmega2[1], logOmega2[nb0.5], logOmega2[nb1], logOmega2[nb4], logOmega2[nb8])


########### ########### ########### ###########

res1 <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res1 <- c(res1, mclapply(1:nbSimu, FUN = one.simu.2stages,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 50,
                             jumpSize = 10,
                             nbK = nbK,
                             varType = "MAD",
                             mc.cores = cores)) ## mc.cores = 8
  }
}
df_2stage <- do.call(rbind, res1)
save(df_2stage, file="df_2stage.RData")

dfmean_1 <- stats::aggregate(df_2stage, list(rep(1:(nrow(df_2stage)%/%nbSimu+1), each = nbSimu, len = nrow(df_2stage))), base::mean)[-1]
z1 <- matrix(dfmean_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2 <- matrix(dfmean_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3 <- matrix(dfmean_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_1 <- stats::aggregate(df_2stage, list(rep(1:(nrow(df_2stage)%/%nbSimu+1), each = nbSimu, len = nrow(df_2stage))), stats::sd)[-1]
w1 <- matrix(dfsd_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2 <- matrix(dfsd_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3 <- matrix(dfsd_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)


########### ########### ########### ###########
########### ########### ########### ###########
########### ########### ########### ###########
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
      colorTable <- designer.colors(nb-1, c( "blue","white"))
      brks <- seq(min, 0, length.out = nb)
    }
    else
    {
      colorTable <- designer.colors(nb-1, c( "white","red"))
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



###
### Plot results
###

plotSimu <- function(z1, w1, z2, w2, z3, w3,
                     .phi = phi,
                     .logOmega2 = logOmega2,
                     .myscale = myscale,
                     .positions = positions)

{
  par(mfrow=c(3,2))
  #par(mar = c(3,4,2,3))
  par(mar = c(3,4,3,3))###change

  z1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
  w1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)

  #####sd_Eta
  u1 <- colScale(min = min(z1), max = max(z1), epsilon = abs(max(z1))/2, nb = 16)
  image.plot(.phi, .logOmega2, z1, breaks=u1$breaks, col=u1$col, axis.args=list(cex.axis=2),
             main = expression(E((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  v1 <- grayScale(max = max(w1),  nb = 20)
  image.plot(.phi, .logOmega2,w1,breaks=v1$breaks, col=v1$col,axis.args=list(cex.axis=2),
             main = expression(SD((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  #####sd_Nu
  u2 <- colScale(min = min(z2), max = max(z2), epsilon = 0.1, nb = 16)
  image.plot(.phi, .logOmega2,z2,breaks=u2$breaks, col=u2$col ,axis.args=list(cex.axis=2),
             main = expression(E((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  v2 <- grayScale(max = max(w2),  nb = 20)
  image.plot(.phi, .logOmega2,w2,breaks=v2$breaks, col=v2$col,axis.args=list(cex.axis=2),
             main = expression(SD((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  #####phi
  u3 <- colScale(min = min(z3), max = max(z3), epsilon = abs(max(z3)/4), nb = 16)
  image.plot(.phi, .logOmega2,z3, breaks=u3$breaks, col=u3$col, axis.args=list(cex.axis=2),
             main = expression(E (hat(phi) - phi) ),
             xlab = expression( phi ),
             ylab = expression( omega^2 ), mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)

  v3 <- grayScale(max = max(w3),  nb = 20)
  image.plot(.phi, .logOmega2,w3,breaks=v3$breaks, col=v3$col,axis.args=list(cex.axis=2),
             main = expression(SD (hat(phi) - phi)),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),  mgp=c(2,2,0),
             cex.lab = 2, cex.main = 2, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95),cex.axis=2)
  axis(2, labels = .myscale, at = .positions,cex.axis=2)
}


plotSimu(z1_1, w1_1, z2_1, w2_1, z3_1, w3_1)


