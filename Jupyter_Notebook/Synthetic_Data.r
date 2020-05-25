devtools::install_github("ykang/gratis")
require("gratis")
library("gratis")
install.packages("rlang")

# Source Paper: https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=SIMULATION+AND+IDENTIFICATION+OF+THE+FRACTIONAL+BROWNIAN+MOTION%3A+A+BIBLIOGRAPHICAL+AND+COMPARATIVE+STUDY&btnG=
#Fractional Brownian Motion Signal:
BM_simulation<-function(n=2000,H=0.7){
  plotfBm <- 1
  
  H2 <- 2 * H
  matcov <- matrix(0, n - 1, n - 1)
  for(i in (1:(n - 1))) {
    j <- i:(n - 1)
    
    r <- 0.5 * (abs(i)^H2 + abs(j)^H2 - abs(j - i)^H2)
    
    r <- r/n^H2
    matcov[i, j] <- r
    matcov[j, i] <- matcov[i, j]
    
  }
  
  
  L <- chol(matcov)
  Z <- rnorm(n - 1,sd=1)
  
  fBm <- t(L) %*% Z
  fBm <- c(0, fBm) ##
  if(plotfBm == 1) {
    par(mfrow = c(1, 1))
    time <- (0:(n - 1))/n
    Nchar <- as.character(n)
    Nleg <- paste(c("N= ", Nchar), collapse = " ")
    Hchar <- as.character(round(H, 3))
    Hleg <- paste(c(", H=", Hchar), collapse = "")
    NHleg <- paste(c(Nleg, Hleg), collapse = "")
    leg <- paste(c("Path of a fractional Brownian motion
      parameters",NHleg), collapse = " : ")
    plot(time, fBm, type = "l", main = leg)
  }
  return (fBm)
}

set.seed(760)
fm_simulted=BM_simulation(n=2000,H=0.7)

#Synthetic Data simulated with gratis-package:
set.seed(123)
# Link to Github-Repository: https://github.com/ykang/gratis/blob/master/man/generate_ts_with_target.Rd
# Link to Paper: https://arxiv.org/pdf/1903.02787.pdf
#n = Number of Time Series to simulate
#ts.length = Number of Samples to simulate for each Time Series
#freq = Frequency of the Time Series
# seasonal: 0 for non-seasonal data, 
#           1 for single-seasonal data, 
#           and 2 for multiple seasonal data
#features = which feature-function to control
#selected.features = for example for stl_features: trend
#target = the desired values for the features

x <- generate_ts_with_target(n = 1, ts.length = 2000, freq = 1, seasonal = 1,
                             features = c('entropy', 'stl_features'),
                             selected.features = c('entropy', 'trend'),
                             target = c(0.3, 0.5))
plot(x)
