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

nbSimu <- 8
nbPhi <- 20
nbOmega2 <- 40
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

# analysis

#dfmean <- aggregate(df,list(rep(1:(nrow(df)%/%nbSimu+1),each=nbSimu,len=nrow(df))),mean)[-1]
#z1 <- matrix(dfmean$`sdEta2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
#image.plot(phi,omega2,z1)
#z2 <- matrix(dfmean$`sdNu2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
#image.plot(phi,omega2,z2)
#z3 <- matrix(dfmean$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
#image.plot(phi,omega2,z3)

#dfsd <- aggregate(df,list(rep(1:(nrow(df)%/%nbSimu+1),each=nbSimu,len=nrow(df))),sd)[-1]
#w1 <- matrix(dfsd$`sdEta2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
#image.plot(phi,omega2,w1)
#w2 <- matrix(dfsd$`sdNu2Est%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
#image.plot(phi,omega2,w2)
#w3 <- matrix(dfsd$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
#image.plot(phi,omega2,w3)

