rm(list = ls())
library(dplyr)

compute.rhoTInf = function(rho, deltaX, deltaY, beta=1/4) {
  rhoTInf1 = rho  + beta * deltaX * deltaY * 0.5 * rho^2;
  rhoTInf2 = 1 + beta * 0.5 * deltaX^2;
  rhoTInf3 = 1 + beta * 0.5 * deltaY^2;
  rhoTInf = rhoTInf1/sqrt(rhoTInf2 * rhoTInf3);
  
}


get_chol <- function(n_gene, rho){
  
  sig <- matrix(rho, nrow = n_gene, ncol = n_gene)
  diag(sig) <- 1
  
  achol <- t(chol(sig))
  
}

## generate correlated normal data 
correlated_norm <- function(mu1, mu2, n1, n2, sigma_chol){
  
  n_gene <- length(mu1)
  y01 <- matrix(rnorm(n_gene*n1), nrow = n_gene, ncol = n1)
  y02 <- matrix(rnorm(n_gene*n2), nrow = n_gene, ncol = n2)
  # leverage the cholesky decoposition to speed up
  y1 <- sigma_chol %*% y01 + mu1
  y2 <- sigma_chol %*% y02 + mu2
  
  dat <- cbind(y1, y2)
  return(dat)
}


## theoretical value of rhot 
rhot <- function(delta, rho){
  
  m <- length(delta)
  rhot_mat <- matrix(NA, ncol = m, nrow = m)

  for(i in 1:(m-1)){
    for(j in (i+1):m)
    rhot_mat[i, j] <- compute.rhoTInf(rho = rho, deltaX = delta[i], deltaY = delta[j], beta = 1/4)
  }
 rho1 <- as.vector(rhot_mat)
 mean(rho1[!is.na(rho1)])
}


# calculate rhot and rhos using simulated data
rhost <- function(dat1, trt){
  
  runlm <- function(y, x){
    x <- as.factor(x)
    lm1 <- summary(lm(y~x))
    # res1 <- lm1$residuals
    # rse <- lm1$sigma
    return(lm1)
  }
  m <- nrow(dat1)
  resid <- matrix(NA, nrow = m, ncol = ncol(dat1))
  
  delta_est <- rep(NA, m)
  for(i in 1:m){
    temp <- runlm(y = dat1[i, ], x = trt)
    resid[i, ] <- temp$residuals
    delta_est[i] <- temp$coefficients[2, 1]
  }
  
  rhot_mat2 <- matrix(NA, ncol = m, nrow = m)
  rhos_mat <- rhot_mat2
  
  for(i in 1:(m-1)){
    for(j in (i+1):m) {
      rhos <- cor(resid[i, ], resid[j, ])
      deltax <- delta_est[i];
      deltay <- delta_est[j]
      rhot_mat2[i, j] <- compute.rhoTInf(rho = rhos, deltaX = deltax, deltaY = deltay, beta = 1/4)
      rhos_mat[i, j] <- rhos
    }
  }
  rho2 <- as.vector(rhot_mat2)
  rhot <- rho2[!is.na(rho2)]
  rho3 <- as.vector(rhos_mat)
  rhos <- rho3[!is.na(rho3)]
  
  res <- data.frame(rhosbar = mean(rhos), rhotbar = mean(rhot))
  return(res)
}


# calculate theoretical vif
vif_theory <- function(delta, rho){
  
  m <- length(delta)
  # theoretical rhot
  theo_rhot <- rhot(delta = delta, rho = rho)
  
  vif_s_th <- 1 + (m-1)*rho
  vif_t_th <- 1 + (m-1)*theo_rhot
  
  return(data.frame(vif_s_th = vif_s_th, vif_t_th = vif_t_th))
}

# calculate vif from simulation
vif_sim <- function(m, n1, n2, delta, rho){
  
  trt <- rep(c(0, 1), c(n1, n2))
  
  sigma_chol <- get_chol(n_gene = m, rho = rho)
  
  dat1 <- correlated_norm(mu1 = rep(0, length(delta)), mu2 = delta, 
                          n1 = n1, n2 = n2, sigma_chol = sigma_chol)
  

  # rhot and rhos from sample 
  rhot_and_rhos <- rhost(dat1 = dat1, trt = trt)
  
  # calculate the vif
  vif_s_sim <- 1 + (m-1) * rhot_and_rhos$rhosbar # using data-row correlation
  vif_t_sim <- 1 + (m-1) * rhot_and_rhos$rhotbar # using calculated rhot correlation
  
  return(data.frame(vif_s_sim = vif_s_sim, vif_t_sim = vif_t_sim))
  
}

#######  perform the simulations

nsim <- 5000 
m <- 21; n1 <- n2 <- 30;
delta = seq(-3, 3, length = m);
rho = 0.1;

vif1 <- vif_theory(delta = delta, rho = rho)
vif2 <- list()
for(i in 1:nsim){
  vif2[[i]] <- vif_sim(m = m, n1 = n1, n2 = n2, delta = delta, rho = rho)
}

vif2 <- bind_rows(vif2)  


## 
# vif_s_sim: vif calculated by simulation using data row correlation
# vif_t_sim: vif calculated by simulation using rho_t estimated from data row correlation
# vif_t_th: vif calculated by formula using rho_t that are related to delta.
# vif_s_th: vif calculated by formula using rho = 0.1
vif3 <- vif2 %>%  mutate(vif_t_th = vif1$vif_t_th, vif_s_th = vif1$vif_s_th)
apply(vif3, 2, mean)
apply(vif3, 2, median)
library(tidyr)
library(ggplot2)
dat <- vif2 %>% gather(key = "source", value = "VIF") %>% 
  mutate(source = ifelse(source == "vif_s_sim", "estimated by data-row correlation", 
                         "estimated by test-statistic correlation"))
ggplot(data = dat, aes(x = VIF, color = source)) + geom_density() + 
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = vif1$vif_t_th, linetype = "dotted", color = "black")

