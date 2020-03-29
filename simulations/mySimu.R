rm(list=ls())
library(ARRWestim)
library(parallel)
cores <- detectCores()
cores <- 45

#########################################################################
#nbSimu <- 10
#lres_1 <- mclapply(1:nbSimu, FUN = one.simu,
#                   n = 10^5,
#                   sdEta2 = 1,
#                   sdNu2 = 1,
#                   phi = 0.5,
#                   poisParam = 0,
#                   meanGap = 2,
#                   nbK = 10,
#                   mc.cores = cores)
#df <- do.call(rbind, lres_1)
#save(df, file="df.RData")
##########################################################################

nbSimu <- 50
nbPhi <- 20
nbOmega2 <- 20
phi <- seq(from = 0, to = 1, length.out = nbPhi+1)
phi <- phi[-(nbPhi+1)]
omega2 <- seq(from = 0.5, to = 10, length.out = nbOmega2)
res <- NULL
for(i in phi)
{
  for(j in omega2)
  {
    res <- c(res, mclapply(1:nbSimu, FUN = one.simu,
                       n = 10^5,
                       sdEta2 = j,
                       sdNu2 = 1,
                       phi = i,
                       poisParam = 0,
                       meanGap = 2,
                       nbK = 10,
                       mc.cores = cores))
    print(i)
    print(j)
  }
}


df <- do.call(rbind, res)
save(df, file="df_sdEta2sdNu2Phi.RData")

# analysis = some simple plots
library(fields)
dfmean <- aggregate(df,list(rep(1:(nrow(df)%/%nbSimu+1),each=nbSimu,len=nrow(df))),mean)[-1]
z1 <- matrix(dfmean$`sdEta2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
image.plot(phi,omega2,z1)
z2 <- matrix(dfmean$`sdNu2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
image.plot(phi,omega2,z2)
z3 <- matrix(dfmean$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
image.plot(phi,omega2,z3)

dfsd <- aggregate(df,list(rep(1:(nrow(df)%/%nbSimu+1),each=nbSimu,len=nrow(df))),sd)[-1]
w1 <- matrix(dfsd$`sdEta2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
image.plot(phi,omega2,w1)
w2 <- matrix(dfsd$`sdNu2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
image.plot(phi,omega2,w2)
w3 <- matrix(dfsd$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
image.plot(phi,omega2,w3)


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

  return(list(col=colorTable, breaks = brks))
}

####################
# generating a 3x2 multi plot with correct color scales

library(fields)
par(mfrow=c(3,2))
par(mar = c(2,2,2,2))

u1 <- colScale(min = min(z1), max = max(z1), epsilon = abs(max(z1))/2, nb = 20)
image.plot(phi,omega2,z1,breaks=u1$breaks, col=u1$col ,axis.args=list(cex.axis=1.3),
           main = expression(E((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])))
v1 <- colScale(min = 0, max = max(w1),  nb = 64)
image.plot(phi,omega2,w1,breaks=v1$breaks, col=v1$col,axis.args=list(cex.axis=1.3),
           main = expression(SD((hat(sigma)[eta]-sigma[eta] )/ sigma[eta])))

u2 <- colScale(min = min(z2), max = max(z2), epsilon = abs(min(z2)/2), nb = 20)
image.plot(phi,omega2,z2,breaks=u2$breaks, col=u2$col ,axis.args=list(cex.axis=1.3),
           main = expression(E((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])))
v2 <- colScale(min = 0, max = max(w2),  nb = 64)
image.plot(phi,omega2,w2,breaks=v2$breaks, col=v2$col,axis.args=list(cex.axis=1.3),
           main = expression(SD((hat(sigma)[nu]-sigma[nu] )/ sigma[nu])))

u3 <- colScale(min = min(z3), max = max(z3), epsilon = abs(max(z3)/2), nb = 20)
image.plot(phi,omega2,z3,breaks=u3$breaks, col=u3$col ,axis.args=list(cex.axis=1.3),
           main = expression(E (hat(phi))))
v3 <- colScale(min = 0, max = max(w3),  nb = 64)
image.plot(phi,omega2,w3,breaks=v3$breaks, col=v3$col,axis.args=list(cex.axis=1.3),
           main = expression(SD (hat(phi))))
