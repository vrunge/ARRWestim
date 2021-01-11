
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
cores <- 1 ###CHANGE

###
### Simulations parameters
###
#nbSimu <- 1200
nbSimu <- 100 ###CHANGE
nbPhi <- 19 #step size = 0.05 in phi
nbOmega2 <- 40

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
  if(min == 0) #variance case
  {
    colorTable <- designer.colors(nb-1, c("white", "black"))
    brks<- seq(0, max, length.out = nb)
  }

  return(list(col = colorTable, breaks = brks))
}

########### ########### ########### ###########

###
### Simulations
###

res <- NULL
for(i in phi)
{
  for(j in omega2)
  {
    res <- c(res, mclapply(1:nbSimu, FUN = one.simu,
                           N = 5000,
                           sdEta = sqrt(j),
                           sdNu = 1,
                           phi = i,
                           nbSeg = 1,
                           jumpSize = 0,
                           nbK = 10,
                           mc.cores = 1)) ## mc.cores = 8
    print(i)
    print(j)
  }
}


###
### Save result
###

df <- do.call(rbind, res)
save(df, file="df_sdEtasdNuPhi_NOCHANGE.RData")

###
### Plot results
###

dfmean <- aggregate(df,list(rep(1:(nrow(df)%/%nbSimu+1),each = nbSimu,len = nrow(df))),mean)[-1]

z1 <- matrix(dfmean$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2 <- matrix(dfmean$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3 <- matrix(dfmean$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd <- aggregate(df,list(rep(1:(nrow(df)%/%nbSimu+1),each = nbSimu,len = nrow(df))),sd)[-1]

w1 <- matrix(dfsd$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2 <- matrix(dfsd$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3 <- matrix(dfsd$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)


par(mfrow=c(3,2))
#par(mar = c(3,4,2,3))
par(mar = c(3,4,2,2))###change


z1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
w1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)

#####sd_Eta
u1 <- colScale(min = min(z1), max = max(z1), epsilon = abs(max(z1))/2, nb = 16)
image.plot(phi, logOmega2, z1, breaks=u1$breaks, col=u1$col, axis.args=list(cex.axis=1.2),
           main = expression(E((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
           xlab = expression( phi ),
           ylab = expression( omega^2 ),   mgp=c(2,2,0),
           cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
           smallplot = c(0.85,0.89,0.13,0.91))
axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
axis(2, labels = myscale, at = positions)

v1 <- colScale(min = 0, max = max(w1),  nb = 20)
image.plot(phi,logOmega2,w1,breaks=v1$breaks, col=v1$col,axis.args=list(cex.axis=1.2),
           main = expression(SD((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])),
           xlab = expression( phi ),
           ylab = expression( omega^2 ),   mgp=c(2,2,0),
           cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
           smallplot = c(0.85,0.89,0.13,0.91))
axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
axis(2, labels = myscale, at = positions)

#####sd_Nu
u2 <- colScale(min = min(z2), max = max(z2), epsilon = abs(min(z2)), nb = 16)
image.plot(phi,logOmega2,z2,breaks=u2$breaks, col=u2$col ,axis.args=list(cex.axis=1.2),
           main = expression(E((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
           xlab = expression( phi ),
           ylab = expression( omega^2 ),   mgp=c(2,2,0),
           cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
           smallplot = c(0.85,0.89,0.13,0.91))
axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
axis(2, labels = myscale, at = positions)

v2 <- colScale(min = 0, max = max(w2),  nb = 20)
image.plot(phi,logOmega2,w2,breaks=v2$breaks, col=v2$col,axis.args=list(cex.axis=1.2),
           main = expression(SD((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])),
           xlab = expression( phi ),
           ylab = expression( omega^2 ),   mgp=c(2,2,0),
           cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
           smallplot = c(0.85,0.89,0.13,0.91))
axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
axis(2, labels = myscale, at = positions)

#####phi
u3 <- colScale(min = min(z3), max = max(z3), epsilon = abs(max(z3)/4), nb = 16)
image.plot(phi, logOmega2,z3, breaks=u3$breaks, col=u3$col, axis.args=list(cex.axis=1.2),
           main = expression(E (hat(phi) - phi) ),
           xlab = expression( phi ),
           ylab = expression( omega^2 ), mgp=c(2,2,0),
           cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
           smallplot = c(0.85,0.89,0.13,0.91))
axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
axis(2, labels = myscale, at = positions)

v3 <- colScale(min = 0, max = max(w3),  nb = 12)
image.plot(phi,logOmega2,w3,breaks=v3$breaks, col=v3$col,axis.args=list(cex.axis=1),
           main = expression(SD (hat(phi) - phi)),
           xlab = expression( phi ),
           ylab = expression( omega^2 ),  mgp=c(2,2,0),
           cex.lab = 1.5, cex.main = 1.5, axes = FALSE,
           smallplot = c(0.85,0.89,0.13,0.91))
axis(1, labels = c(0,0.2,0.4,0.6,0.8,0.95), at = c(0,0.2,0.4,0.6,0.8,0.95))
axis(2, labels = myscale, at = positions)

