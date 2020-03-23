
#library(parallel)
#cores <- detectCores()
#cores <- 40

#############
############# one.simu
#############
one.simu <- function(i, n = 10^5, phi = 0.5, sdEta2 = 1, sdNu2 = 1, nbK, poisParam = 0, meanGap = 2)
{
  y <- dataARRW(n = n, poisParam, meanGap, phi = phi, sdEta2 = sdEta2, sdNu2 = sdNu2)
  v <- estimVar(y = y$y, nbK = nbK) ##vector of estimated variances
  res <- bestParameters(v = v, nbK = nbK)

  #ind, method, mse, beta, K
  df <- data.frame(numeric(0), numeric(0), numeric(0), numeric(0) ,stringsAsFactors = FALSE)
  colnames(df) <- c("index", "sdEta2", "sdNu2", "phi")
  df[1,] <- c(i,(res$Eta2Opt - sdEta2)/sdEta2, (res$Nu2Opt-sdNu2)/sdNu2, res$argmin - phi)
  return(df)
}


#########################################################################
#nbSimu <- 100
#df_test <- NULL
#for(i in 1:10)
#{
#  lres_1 <- mclapply(1:nbSimu, FUN = one.simu,
#                      n = 10^5,
#                      phi = 0.1,
#                     sdEta2 = 1,
#                      sdNu2 = i,
#                      nbK = 10,
#                      mc.cores = cores)
#df_test <- rbind(df_test, do.call(rbind, lres_1))
#}
#df_test



#save(df_test, file="df_test.Rdata")
