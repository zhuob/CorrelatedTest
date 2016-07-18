### simulation 

library(MASS)
library(reshape2)
library(ggplot2)
library(directlabels)  # add labels to the contour plot


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
##  arguments are objects returned from compare_rho function

comb1 <- function(obj1, obj2, obj3){
  
  rho_mat1 <- obj1$rho_mat;
  rho_mat2 <- obj2$rho_mat
  rho_mat1$s2 <- rho_mat2$test_cor
  rho_mat1$theoretical <- obj3
  
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



plot_rho1 <- function(obj, case, textsize = rep(20, 4),
                      smaller =5, larger=1000, include_invS = F){
  
  rho_t <- obj$rho_mat
  rho_s <- obj$rho_sample
  DE <- obj$DE
  
#  rho_t$rho_sample_ave <- apply(rho_s, 1, mean) 
#  names(rho_t)[5] <- "sample"
  v1 <- paste("n =", smaller)
  v2 <- paste("n =", larger)
  names(rho_t)[c(2, 4)] <- c(v1, v2)
  
  if (include_invS){
    prep_data <- rho_t
  } else {
    prep_data <- rho_t[, c(1, 2, 4, 5)]
  }
  
  prep_data2 <- melt(prep_data, id = "true_popu_rho", 
                     variable.name = "category", value.name = "estimate")
  
  tl <- substitute(case1~"):"~"("~delta[x]~","~delta[y]~")"~"="~"("~var1~","~var2~")",
                   list(case1 = case, var1 = DE[1], var2=DE[2]))
  #tl <- substitute(Delta[x]~ "="~ var1 ~","~Delta[y]~"="~var2, list (var1 = DE[1], var2=DE[2]))
  p1 <- ggplot(data = prep_data2, aes(true_popu_rho, estimate, color = category)) + 
    geom_line(aes(linetype = category)) + 
    labs( x= "population correlation", y = "test statistics correlation", 
          title = tl ) + 
    theme(legend.position = c(0.7, 0.2), 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]),
          axis.text = element_text(size = textsize[3]), 
          axis.title = element_text(size = textsize[4], face= "bold"))
  p1  + guides(
    linetype  = guide_legend(keywidth = 3, keyheight = 1), 
    color = guide_legend(keywidth = 3, keyheight = 1) ) +
    geom_abline(intercept = 0, slope=1) + ylim(-1, 1)
    
  
}


## theoretical correlation 
## delta1 = Delta_x/sigma_x;  delta2 = Delta_y/sigma_y

theore_rho <- function(delta1, delta2, rho, sigma, n1, n2){
  
  beta <- 1/(4 + 2*n1/n2 + 2*n2/n1)
  delta <- c(delta1, delta2)
  signal2noise <- delta/sigma 
  
  numerator <- rho + beta*prod(signal2noise)*rho^2
  denominator <- sqrt((1 + beta*signal2noise[1]^2)*(1 + beta*signal2noise[2]^2))
  
  rho_stat <- numerator/denominator
  return(rho_stat)
}


theoretical_plot <- function(delta1_vector, delta2= seq(0,4), rho, sigma, n1, n2, textsize = rep(20, 4)){
  
  m <- length(delta2)
  obj_dat <- data.frame(matrix(NA, length(delta1_vector), m))
  for ( k in 1:m)
  {
   obj_dat[, k] <- sapply(delta1_vector, theore_rho, rho = rho, delta2=delta2[k], sigma= sigma, n1=n1, n2=n2)
  }
  names(obj_dat) <- delta2
  obj_dat$delta1 <- delta1_vector
  
  
  tr_data <- melt(obj_dat, id = "delta1", value.name = "theoretical",
                  variable.name = "delta2" )
  rho1 <- rho
  tl <- substitute(rho~"="~rho1, list( rho1 = rho1))
  
  ggplot(tr_data, aes(delta1, theoretical, color = delta2, linetype = delta2)) + 
    geom_line(size = 1) + 
    labs(y = "theoretical correlation", title = tl) + 
    theme(legend.position = c(0.9, 0.2), 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text = element_text(size = textsize[3]), 
          axis.title = element_text(size = textsize[4], face = "bold")) + 
    guides(fill = guide_legend(keywidth = 1, keyheight = 1), 
        #   linetype = guide_legend(keywidth = 3, keywidth= 1), 
           colour = guide_legend(keywidth = 3, keyheight = 1))
    

}


contour_plot <- function(delta1=seq(-4, 4, 0.1), delta2= seq(-4, 4, 0.1), rho, 
                         sigma, n1, n2, brks, textsize=rep(20, 4)){
  
  m1 <- length(delta1)
  m2 <- length(delta2)
  obj_data <- data.frame(matrix(NA, nrow = m1*m2, ncol = 3))
  for ( k in 1: m2){
    start_index <- (k-1)*m1 + 1
    end_index <- k*m1
    obj_data[start_index:end_index, 3] <- sapply(delta1, theore_rho, rho = rho, delta2=delta2[k], sigma= sigma, n1=n1, n2=n2)
    obj_data[start_index:end_index, 1] <- delta1
    obj_data[start_index:end_index, 2] <- delta2[k]
  }
  rho1 <- rho
  tl <- substitute(rho~"="~rho1, list( rho1 = rho1)) # the title

  names(obj_data) <- c("delta.x", "delta.y","TestCorr")
  v <- ggplot(obj_data, aes(delta.x, delta.y, z =TestCorr )) +
    geom_contour(aes(colour=..level..), breaks = brks) + 
    labs(title = tl, x = expression(delta[x]), y = expression(delta[y])) + 
    theme(legend.position = "none",
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text = element_text(size = textsize[3]),
          axis.title = element_text(size = textsize[4], face = "bold")
          )
    direct.label(v,method="bottom.pieces")

}


