rm(list=ls())



library(ARRWestim)
library(parallel)
cores <- detectCores()
cores <- 1

#########################################################################
nbSimu <- 10
lres_1 <- mclapply(1:nbSimu, FUN = one.simu,
                   n = 10^5,
                   sdEta2 = 1,
                   sdNu2 = 1,
                   phi = 0.5,
                   poisParam = 0,
                   meanGap = 2,
                   nbK = 10,
                   mc.cores = cores)

df_linear_gauss <- do.call(rbind, lres_1)

#save(df_linear_gauss, file="df_linear_gauss.RData")
