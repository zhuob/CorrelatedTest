### simulation 

library(MASS)
library(reshape2)
library(ggplot2)

correlated.norm <- function(mu, rho, sigma, n){

## generate correlated Normal data
## input:  
#    mu:  a vector of 2, specifying the mean level for each gene
#   rho:  the true correlation we want to simulate 
# sigma:  a vector of 2, specifying the variance 
#     n:  the number of samples to be generated

## output: a matrix of 2 by n

  sig.element <- c(sigma[1], rho*sqrt(sigma[1]*sigma[2]),     ## the elements in Sigma
                   rho*sqrt(sigma[1]*sigma[2]), sigma[2])
  sig <- matrix(sig.element, 2, 2)                            # the covariance matrix
  dat <- mvrnorm(n=n, mu = mu, Sigma = sig)                   # simulate the correlated normal data
  return(t(dat))
}



ttest.stat <- function(y){
  
# calculate the t-test statistics 
# input:  
#  y: 2 by n matrix
# output:
#  a vector of two t-statistics 
  
  n <- dim(y)[2]                              # get the number of samples from each gene
  mean1 <- apply(y[, 1:(n/2)], 1, mean)       # the first half is from treatment
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)    # the second half is from control
  
  s1 <- apply(y[, 1:(n/2)], 1, var)           # the variance for the first half
  s2 <- apply(y[, -(1:(n/2))], 1, var)        # the variance for the second half

#   t1 <- t.test(y[1, 1:(n/2)], y[1, -(1:(n/2))])$stat
#   t2 <- t.test(y[2, 1:(n/2)], y[2, -(1:(n/2))])$stat
#   t.stat <- c(t1, t2)

  t.stat <- (mean1 - mean2)/sqrt(s1/(n/2) + s2/(n/2))  # the test statistics
  return(t.stat)
}



ztest.stat <- function(y, sigma1, sigma2){
  
# calculate the z-test statistics 
# input:  
#  y: 2 by n matrix
# output:
#  a vector of two t-statistics 
  
  n <- dim(y)[2]                              # get the number of samples from each gene
  mean1 <- apply(y[, 1:(n/2)], 1, mean)       # the first half is from treatment
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)    # the second half is from control
  
  z.stat <- (mean1-mean2)/sqrt(sigma1/(n/2) + sigma2/(n/2))
  return(z.stat)
}




# simulate the correlation of 1/s1 and 1/s2

##  calculate the inverse of standard errors
inv_sd <- function(y){
  
  n <- dim(y)[2]                              # get the number of samples from each gene
  var1 <- apply(y[, 1:(n/2)], 1, var)         # the variance for the first half
  var2 <- apply(y[, -(1:(n/2))], 1, var)      # the variance for the second half
  
  inv_s1 <- 1/sqrt((var1[1] + var2[1])/2)                 # pooled variance for gene 1
  inv_s2 <- 1/sqrt((var1[2] + var2[2])/2)                 # pooled variance for gene 2
  
  #  inv_s1 <- (var1[1] + var2[1])/2                 # pooled variance for gene 1
  #  inv_s2 <- (var1[2] + var2[2])/2                 # pooled variance for gene 2
  
  return(c(inv_s1, inv_s2))
}





############ now the correlation between test statistic is related to both 
###########  rho0 = cor(z1, z2) and rho1 = cor(1/s1, 1/s2)

# the correlation of sample and the correlation of inverse_standard_errors
sample.inv_sd.correlation <- function(mu1, sigma1, mu2, sigma2, rho, n, nreps){

  stat <- matrix(NA, nrow=nreps, ncol=2)                     # matrix to store the t.stat
  rho.sample <- c()
  
  for ( k in 1: nreps)                                       # create a vector of score stat for each gene
  {
    y.t1 <- correlated.norm(mu1, rho, sigma1, n/2)           # treatment 
    y.t2 <- correlated.norm(mu2, rho, sigma2, n/2)           # control
    
    y <- cbind(y.t1, y.t2)                                   # expression data of two genes
    rho1 <- cor(t(y.t1))[1, 2]; rho2 <- cor(t(y.t2))[1, 2]   # sample correlation for each group
    rho.sample[k] <- mean(c(rho1, rho2))                     # mean of sample correlation
    
    #    stat[k, ] <- ztest.stat(y, sigma1, sigma2)              # calculate the z stat
    stat[k, ] <- inv_sd(y)
  }
  
  rho.ave <- mean(rho.sample)                                # calculate the mean of the sampel rhos
  rho.s <- cor(stat)[1, 2]                                   # calculate the correlation of statistics
  rho.list <- list(sample.cor= rho.ave, inv_sd.cor = rho.s)
  
  return(rho.list)  
}





sample.stat.correlation <- function(mu1, sigma1, mu2, sigma2, rho, n, nreps, test.type = "t"){
#   simulate correlated normal data as specified by the parameters and 
#  return the sample correlation as  well as the test stat correlation

#  input:
#   mu1: true mean for group 1
#   mu2: true mean for group 2
# sigma: true variance 
#   rho: true correlation
#     n: total number of samples to be simulated
# nreps: how many data sets to be simulated
#
# output:
#  rho.list: the estimated sample correlation and estimated test statistics correlation
  
  stat <- matrix(NA, nrow=nreps, ncol=2)                     # matrix to store the t.stat
  rho.sample <- c()
  
  for ( k in 1: nreps)                                       # create a vector of score stat for each gene
  {
    y.t1 <- correlated.norm(mu1, rho, sigma1, n/2)           # treatment 
    y.t2 <- correlated.norm(mu2, rho, sigma2, n/2)           # control
    
    y <- cbind(y.t1, y.t2)                                   # expression data of two genes
    rho1 <- cor(t(y.t1))[1, 2]; rho2 <- cor(t(y.t2))[1, 2]   # sample correlation for each group
    rho.sample[k] <- mean(c(rho1, rho2))                     # mean of sample correlation
    
    if (test.type == "t") { 
      stat[k, ] <- ttest.stat(y)                               #  calculate the t stat
      }
    else if ( test.type == "z"){
      stat[k, ] <- ztest.stat(y, sigma1, sigma2)              # calculate the z stat
    }
    
  }
  
  rho.ave <- mean(rho.sample)                                # calculate the mean of the sampel rhos
  rho.t <- cor(stat)[1, 2]                                   # calculate the correlation of statistics
  rho.list <- list(true.cor = rho, sample.cor= rho.ave, stat.cor = rho.t)
  
  return(rho.list)  
}


simulate.rho <- function(mu1, sigma1, mu2, sigma2, rhoVec, n, nreps, test.type= "t"){
# simulate sample and test statistics correlation according to the rho Vector 
  
  rho <- rhoVec
  store.corr <- data.frame(matrix(NA, length(rho), 3))
  
  for ( i in 1:length(rho))
  {
    a <- sample.stat.correlation(mu1, sigma1, mu2, sigma2, rho[i], n, nreps, test.type = test.type)
    store.corr[i, ] <- unlist(a)  
  }
  colnames(store.corr) <- c("true", "sample", "test.stat")
  
  return(store.corr)
}


plot.rho <- function(data, textsize=rep(20, 4)){
# plot the sample and test statistics correlation against the true correlation  

  store.corr1 <- melt(data, id="true", value.name = "correlation", variable.name = "type")

    ggplot(data = store.corr1, aes(x=true, y=correlation) ) +
      geom_point(aes(colour=type, shape=type), size=2) + 
      #labs(x="true", y="correlation") + 
      theme(legend.position="top", 
            legend.text = element_text(size = textsize[1]),
            plot.title = element_text(size = textsize[2]), 
            axis.text=element_text(size=textsize[3]), 
            axis.title=element_text(size=textsize[4],face="bold"))
  
}


