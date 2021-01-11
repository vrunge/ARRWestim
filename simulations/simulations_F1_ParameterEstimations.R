
rm(list=ls())

# LIBRARIES
library(ARRWestim)
library(parallel)
library(fields)

###
### NB cores for parallel computing
###
cores <- detectCores()
#cores <- 45
cores <- 8 ###CHANGE

###
### Simulations parameters
###
#nbSimu <- 1200
nbSimu <- 1000 ###CHANGE
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
  if(min*max < 0) #biais estimation case
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
  if(min == 0) #variance case
  {
    colorTable <- designer.colors(nb-1, c("white", "black"))
    brks<- seq(0, max, length.out = nb)
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
  for(j in omega2)
  {
    res1 <- c(res1, mclapply(1:nbSimu, FUN = one.simu,
                           N = 5000,
                           sdEta = sqrt(j),
                           sdNu = 1,
                           phi = i,
                           nbSeg = 1,
                           jumpSize = 0,
                           nbK = nbK,
                           mc.cores = 8)) ## mc.cores = 8
  }
  print(i)
}





res2 <- NULL
for(i in phi)
{
  for(j in omega2)
  {
    res2 <- c(res2, mclapply(1:nbSimu, FUN = one.simu,
                           N = 5000,
                           sdEta = sqrt(j),
                           sdNu = 1,
                           phi = i,
                           type = "rand1",
                           nbSeg = 50,
                           jumpSize = 10,
                           nbK = nbK,
                           mc.cores = 8)) ## mc.cores = 8
  }
  print(i)
}


res3 <- NULL
for(i in phi)
{
  for(j in omega2)
  {
    res3 <- c(res3, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 100,
                             jumpSize = 10,
                             nbK = nbK,
                             mc.cores = 8)) ## mc.cores = 8
  }
  print(i)
}
###
### Save result
###


df <- do.call(rbind, res1)
save(df, file="df_sdEtasdNuPhi_noCHANGE.RData")

dfr1 <- do.call(rbind, res2)
save(dfr1, file="df_sdEtasdNuPhi_rand1_50CHANGES.RData")

dfr2 <- do.call(rbind, res3)
save(dfr2, file="df_sdEtasdNuPhi_rand1_100CHANGES.RData")


dfmean_1 <- stats::aggregate(df, list(rep(1:(nrow(df)%/%nbSimu+1), each = nbSimu, len = nrow(df))), base::mean)[-1]
z1_1 <- matrix(dfmean_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2_1 <- matrix(dfmean_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3_1 <- matrix(dfmean_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_1 <- stats::aggregate(df, list(rep(1:(nrow(df)%/%nbSimu+1), each = nbSimu, len = nrow(df))), stats::sd)[-1]
w1_1 <- matrix(dfsd_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2_1 <- matrix(dfsd_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3_1 <- matrix(dfsd_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)


dfmean_2 <- stats::aggregate(dfr1, list(rep(1:(nrow(dfr1)%/%nbSimu+1), each = nbSimu, len = nrow(dfr1))), base::mean)[-1]
z1_2 <- matrix(dfmean_2$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2_2 <- matrix(dfmean_2$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3_2 <- matrix(dfmean_2$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_2 <- stats::aggregate(dfr1, list(rep(1:(nrow(dfr1)%/%nbSimu+1), each = nbSimu, len = nrow(dfr1))), stats::sd)[-1]
w1_2 <- matrix(dfsd_2$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2_2 <- matrix(dfsd_2$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3_2 <- matrix(dfsd_2$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)


dfmean_3 <- stats::aggregate(dfr2, list(rep(1:(nrow(dfr2)%/%nbSimu+1), each = nbSimu, len = nrow(dfr2))), base::mean)[-1]
z1_3 <- matrix(dfmean_3$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2_3 <- matrix(dfmean_3$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3_3 <- matrix(dfmean_3$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_3 <- stats::aggregate(dfr2, list(rep(1:(nrow(dfr2)%/%nbSimu+1), each = nbSimu, len = nrow(dfr2))), stats::sd)[-1]
w1_3 <- matrix(dfsd_3$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2_3 <- matrix(dfsd_3$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3_3 <- matrix(dfsd_3$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)


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
  par(mar = c(3,4,2,2))###change

  z1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
  w1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)

  #####sd_Eta
  u1 <- colScale(min = min(z1), max = max(z1), epsilon = abs(max(z1))/2, nb = 16)
  image.plot(.phi, .logOmega2, z1, breaks=u1$breaks, col=u1$col, axis.args=list(cex.axis=1.2),
             main = expression(E((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
  axis(2, labels = .myscale, at = .positions)

  v1 <- grayScale(max = max(w1),  nb = 20)
  image.plot(.phi, .logOmega2,w1,breaks=v1$breaks, col=v1$col,axis.args=list(cex.axis=1.2),
             main = expression(SD((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
  axis(2, labels = .myscale, at = .positions)

  #####sd_Nu
  u2 <- colScale(min = min(z2), max = max(z2), epsilon = 0.1, nb = 16)
  image.plot(.phi, .logOmega2,z2,breaks=u2$breaks, col=u2$col ,axis.args=list(cex.axis=1.2),
             main = expression(E((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
  axis(2, labels = .myscale, at = .positions)

  v2 <- grayScale(max = max(w2),  nb = 20)
  image.plot(.phi, .logOmega2,w2,breaks=v2$breaks, col=v2$col,axis.args=list(cex.axis=1.2),
             main = expression(SD((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),   mgp=c(2,2,0),
             cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
  axis(2, labels = .myscale, at = .positions)

  #####phi
  u3 <- colScale(min = min(z3), max = max(z3), epsilon = abs(max(z3)/4), nb = 16)
  image.plot(.phi, .logOmega2,z3, breaks=u3$breaks, col=u3$col, axis.args=list(cex.axis=1.2),
             main = expression(E (hat(phi) - phi) ),
             xlab = expression( phi ),
             ylab = expression( omega^2 ), mgp=c(2,2,0),
             cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
  axis(2, labels = .myscale, at = .positions)

  v3 <- grayScale(max = max(w3),  nb = 20)
  image.plot(.phi, .logOmega2,w3,breaks=v3$breaks, col=v3$col,axis.args=list(cex.axis=1),
             main = expression(SD (hat(phi) - phi)),
             xlab = expression( phi ),
             ylab = expression( omega^2 ),  mgp=c(2,2,0),
             cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
             smallplot = c(0.85,0.89,0.13,0.91))
  axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
  axis(2, labels = .myscale, at = .positions)
}


plotSimu(z1_1, w1_1, z2_1, w2_1, z3_1, w3_1)
plotSimu(z1_2, w1_2, z2_2, w2_2, z3_2, w3_2)
plotSimu(z1_3, w1_3, z2_3, w2_3, z3_3, w3_3)
