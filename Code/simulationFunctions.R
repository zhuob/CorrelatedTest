### simulation 

library(MASS)
library(reshape2)
library(ggplot2)


#' generate bivariate normally distributed data.
#'
#' @title simulate one gene expression experiment for two group comparison 
#' @param mu1  a vector of 2, specifying the means in Group 1
#' @param  rho  the true correlation we want to simulate 
#' @param sigma  a vector of 2, specifying the population variance 
#' @param n1  the number of samples to be generated in group 1
#' @param n2 the number of samples in group 2
#' @return a list 
#' \item{y}{a matrix of 2 by n, where n = n1 + n2}
#' \item{trt}{treatment labels, 1 for group 1 and 2 for group 2}
#' @export

correlated_norm <- function(mu1, mu2, rho, sigma, n1, n2){
  sig_element <- c(sigma[1], rho*sqrt(sigma[1]*sigma[2]),     ## the elements in Sigma
                   rho*sqrt(sigma[1]*sigma[2]), sigma[2])
  sig <- matrix(sig_element, 2, 2)                            # the covariance matrix
  y1 <- mvrnorm(n=n1, mu = mu1, Sigma = sig)                   # simulate the correlated normal data for group 1
  y2 <- mvrnorm(n=n2, mu = mu2, Sigma = sig)
  y <- rbind(y1, y2)
  trt <- c(rep(1, n1), rep(2, n2))
  dat <- list(data = t(y), trt = trt)
  return(dat)
}


#' calculate t test statistics from the data.
#'
#' @title t test statistics for a pair of simulated genes 
#' @param dat  an object returned from correlated_norm function
#' @return a list
#' \item{t_stat}{a vector of two, the t-test statistics}
#' \item{r}{estimated sample correlation}
#' \item{spool}{a vector of two, the corresponding pooled variances}
#' @export

ttest_stat <- function(dat){
  
# calculate the t-test statistics 
# input:  
#  y: 2 by n matrix
# output:
#  a vector of two t-statistics 
  y <- dat$data
  trt <- dat$trt
  y1 <- y[, trt==1]
  y2 <- y[, trt!=1]
  n1 <- ncol(y1)  
  n2 <- ncol(y2)
  
  mean1 <- apply(y1, 1, mean)
  mean2 <- apply(y2, 1, mean)
  s1 <- apply(y1, 1, var)
  s2 <- apply(y2, 1, var)
  
  sp <- ((n1-1)*s1 + (n2-1)*s2) / (n1 + n2-2)
  

#   t1 <- t.test(y[1, 1:(n/2)], y[1, -(1:(n/2))])$stat
#   t2 <- t.test(y[2, 1:(n/2)], y[2, -(1:(n/2))])$stat
#   t.stat <- c(t1, t2)
  
  t_stat <- (mean1 - mean2)/sqrt(sp * (1/n1 + 1/n2))  # the test statistics
  
  y_resid <- cbind(y1-mean1, y2-mean2)
  r <- cor(t(y_resid))[1, 2]
  
  results <- list(t_stat= t_stat, r = r, spool = sp)
  
  return(results)
}


#' repeat one simulation to get a multiple pairs of test statistics
#'
#' @title multiple test statistics
#' @param mu1  a vector of 2, specifying the means in Group 1
#' @param  rho  the true correlation we want to simulate 
#' @param sigma  a vector of 2, specifying the population variance 
#' @param n1  the number of samples to be generated in group 1
#' @param n2 the number of samples in group 2
#' @param nrep   number of experiment to be generated
#' @return a data frame 
#' @export

## simulate the test statistics correlation 
test_cor <- function(mu1, mu2, rho, sigma, n1, n2, nreps){
  
  stat <- matrix(NA, nrow= nreps, ncol=2)
  r_sample <- matrix(NA, nrow= nreps, ncol=1)
  inv_s <- stat
  
  for (k in 1: nreps)
  {
    dat <- correlated_norm(mu1, mu2, rho, sigma, n1, n2)
    result <- ttest_stat(dat)
    stat[k, ] <- result$t_stat          
    inv_s[k, ] <- 1/sqrt(result$spool)  # calculate s1^{-1} and s2^{-1}
    r_sample[k, ] <- result$r
  }
  
  results <- data.frame(stat, inv_s, r_sample)
  names(results) <- c("stat1", "stat2", "inv_s1", "inv_s2", "sample_cor")
  return(results)

}




#' simulate test statistics correlation for a serial of rhos
#'
#' @title comparison of correlations
#' @param mu1  a vector of 2, specifying the means in Group 1
#' @param  rhovec  the true correlation vector we want to simulate 
#' @param sigma  a vector of 2, specifying the population variance 
#' @param n1  the number of samples to be generated in group 1
#' @param n2 the number of samples in group 2
#' @param nrep   number of experiment to be generated
#' @return a list
#' \item{rho_mat}{a data frame }
#' \item{rho_sample}{a matrix}

#' @export

compare_rho <- function(mu1, mu2, rhovec, sigma, n1, n2, nreps){
  
  
  rho_mat <- matrix(NA, nrow = length(rhovec), 3)
  rho_sample <- matrix(NA, nrow= length(rhovec), ncol = nreps)
  
  rho_mat <- data.frame(rho_mat)
  names(rho_mat) <- c("true_popu_rho", "test_cor", "invS_cor")
  
  for ( k in 1:length(rhovec))
  {
    single_rho <- test_cor(mu1, mu2, rhovec[k], sigma, n1, n2, nreps)
    a1 <- cor(single_rho[, 1:2])[1, 2]
    a2 <- cor(single_rho[, 3:4])[1, 2]
    rho_mat[k, ] <- c(rhovec[k], a1, a2)
    rho_sample[k, ] <- as.vector(single_rho[, 5])
  }
  
  ensemble <- list(rho_mat = rho_mat, rho_sample=rho_sample, DE = mu2-mu1)
  return(ensemble)
  
}


#' plot the simulated rhos
#'
#' @title plot the simulated rhos
#' @param obj   an object returd by compare_rho function
#' @param  rhovec  the true correlation vector we want to simulate 
#' @param sigma  a vector of 2, specifying the population variance 
#' @param n1  the number of samples to be generated in group 1
#' @param n2 the number of samples in group 2
#' @param nrep   number of experiment to be generated
#' @return a list
#' \item{rho_mat}{a data frame }
#' \item{rho_sample}{a matrix}

#' @export

plot_rho <- function(obj, case, textsize = rep(20, 4),include_invS = F){
  
  rho_t <- obj$rho_mat
  rho_s <- obj$rho_sample
  DE <- obj$DE
  
  rho_t$rho_sample_ave <- apply(rho_s, 1, mean) 
  names(rho_t)[2:4] <- c("TestStat", "InvS", "Sample")
  
  if (include_invS){
    prep_data <- rho_t
  } else {
    prep_data <- rho_t[, c(1, 2, 4)]
  }
  
  prep_data2 <- melt(prep_data, id = "true_popu_rho", 
                     variable.name = "category", value.name = "estimate")
  
  tl <- substitute(case1~"):"~"("~Delta[x]~","~Delta[y]~")"~"="~"("~var1~","~var2~")", 
                   list(case1 = case, var1 = DE[1], var2=DE[2])  )
  #tl <- substitute(Delta[x]~ "="~ var1 ~","~Delta[y]~"="~var2, list (var1 = DE[1], var2=DE[2]))
  p1 <- ggplot(data = prep_data2, aes(true_popu_rho, estimate, color = category)) + 
    geom_line(aes(linetype = category)) + 
    labs( x= "population correlation", y = "estimated correlation", 
          title = tl ) + 
    theme(legend.position = c(0.8, 0.1), 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]),
          axis.text = element_text(size = textsize[3]), 
          axis.title = element_text(size = textsize[4], face= "bold")
          )
  p1  + guides(
      linetype  = guide_legend(keywidth = 3, keyheight = 1), 
      color = guide_legend(keywidth = 3, keyheight = 1)
  )
  
}







##  compare the performance of small samples and large samples
comb1 <- function(obj1, obj2){
  
  rho_mat1 <- obj1$rho_mat;
  rho_mat2 <- obj2$rho_mat
  rho_mat1$s2 <- rho_mat2$test_cor
  
  obj_new <- obj1;
  obj_new$rho_mat <- rho_mat1
  names(obj_new$rho_mat)[2:4] <- c("smaller", "invS", "larger")
  
  return(obj_new)
}


plot_rho2 <- function(obj, case, textsize = rep(20, 4),include_invS = F){
  
  rho_t <- obj$rho_mat
  rho_s <- obj$rho_sample
  DE <- obj$DE
  
  rho_t$rho_sample_ave <- apply(rho_s, 1, mean) 
 names(rho_t)[5] <- "sample"
  
  if (include_invS){
    prep_data <- rho_t
  } else {
    prep_data <- rho_t[, c(1, 2, 4, 5)]
  }
  
  prep_data2 <- melt(prep_data, id = "true_popu_rho", 
                     variable.name = "category", value.name = "estimate")
  
  tl <- substitute(case1~"):"~"("~Delta[x]~","~Delta[y]~")"~"="~"("~var1~","~var2~")", 
                   list(case1 = case, var1 = DE[1], var2=DE[2]))
  #tl <- substitute(Delta[x]~ "="~ var1 ~","~Delta[y]~"="~var2, list (var1 = DE[1], var2=DE[2]))
  p1 <- ggplot(data = prep_data2, aes(true_popu_rho, estimate, color = category)) + 
    geom_line(aes(linetype = category)) + 
    labs( x= "population correlation", y = "estimated correlation", 
          title = tl ) + 
    theme(legend.position = c(0.8, 0.1), 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]),
          axis.text = element_text(size = textsize[3]), 
          axis.title = element_text(size = textsize[4], face= "bold"))
  p1  + guides(
    linetype  = guide_legend(keywidth = 3, keyheight = 1), 
    color = guide_legend(keywidth = 3, keyheight = 1) )
  
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


