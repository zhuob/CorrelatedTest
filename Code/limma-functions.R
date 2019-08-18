# rm(list = ls())
library(dplyr)
library(limma)
library(ggplot2)
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/product.moments.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/plot.rhoT.vs.rho.R")


get_chol <- function(n_gene, rho){
  
  sig <- matrix(rho, nrow = n_gene, ncol = n_gene)
  diag(sig) <- 1
  
  achol <- t(chol(sig))
  
}


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


prepare_mu <- function(n_gene, p_low, p_high, 
                       de_low, de_high){
  
  n_de_low <- n_gene * p_low
  n_de_high <- n_gene * p_high
  n_de <- n_de_low + n_de_high
  
  mu1 <- rep(0, n_gene)
  mu2 <- rep(0, n_gene)
  #mu2[1:n_de_low] <- rbinom(n_de_low, de_low, 0.7) * (-1)
  #mu2[(n_de_low + 1):(n_de_low + n_de_high)] <- rbinom(n_de_high, de_high, 0.8)
  mu2[1:n_de_low] <- de_low
  mu2[(n_de_low + 1):(n_de_low + n_de_high)] <- de_high
  
  return(tibble(mu1=mu1, mu2 = mu2))
}

t_test <- function(dat, trt, row_num){
  
  y1 <- dat[row_num, trt == 0]
  y2 <- dat[row_num, trt == 1]
  t.test(y1, y2, var.equal = TRUE)$statistic
}


simulate_t <- function(nsim, seed = 123, rho, n1, n2, 
                       de_low, de_high, n_gene) {
  
  trt <- c(rep(0, n1), rep(1, n2))
  
  sigma_chol <- get_chol(n_gene = n_gene, rho = rho)
  set.seed(seed)
  mu <- prepare_mu(n_gene, de_low = de_low, de_high = de_high,
                   p_low = 0.05, p_high = 0.05)
  
  sim1 <- data.frame(matrix(NA, nrow = n_gene, ncol = nsim)) # limma statistics
  sim2 <- sim1 # two sample t-test
  names(sim1) <- names(sim2) <- paste("sim", 1:nsim, sep = "_")
  
  time1 <- Sys.time()
  for(i in 1:nsim){
    cat(i, "\r")
    dat1 <- correlated_norm(mu1 = mu$mu1, mu2 = mu$mu2, 
                            n1 = n1, n2 = n2, sigma_chol)
    
    design <- model.matrix(~trt)
    fit <- lmFit(dat1, design = design) %>% ebayes()
    sim1[, i] <- fit$t[, 2]
    t_stat <- purrr::map_dbl(.x = 1:n_gene, t_test, dat = dat1, trt = trt)
    sim2[, i] <- t_stat
  }
  Sys.time() - time1
  
  rslt <- list(t_limma = sim1, tt = sim2)
  return(rslt)
}



plot_t_cor <- function(sim1, rho){
  
  temp <- cor(t(sim1))
  n_gene <- ncol(temp)
  n_de_low <- n_gene*0.05;
  n_de_high <- n_gene*0.05
  n_de <- n_de_low + n_de_high
  
  b1 <- temp[1:n_de_low, 1:n_de_low]  %>% as.vector  # low DE
  b2 <- temp[(n_de_low+1):n_de, (n_de_low+1):n_de] %>% as.vector # high DE
  b3 <- temp[(n_de+1):n_gene, (n_de+1):n_gene] %>% as.vector  # no DE
  # correlation between de genes
  b12 <- temp[1:n_de_low, (n_de_low+1):n_de] %>% as.vector
  b23 <- temp[(n_de_low+1):n_de, (n_de+1):n_gene] %>% as.vector
  b13 <- temp[1:n_de_low, ((n_de +1):n_gene)] %>% as.vector
  
  # remove diagonal elements
  b1 <- b1[b1<1] %>% as.vector %>% as_tibble %>% mutate(cluster = "low DE")
  b2 <- b2[b2<1] %>% as.vector %>% as_tibble %>% mutate(cluster = "high DE")
  b3 <- b3[b3<1] %>% as.vector %>% as_tibble %>% mutate(cluster = "Non-DE")
  b12 <- b12 %>% as.vector %>% as_tibble %>% mutate(cluster = "low vs high DE")
  b23 <- b23 %>% as.vector %>% as_tibble %>% mutate(cluster = "high vs Non DE")
  b13 <- b13 %>% as.vector %>% as_tibble %>% mutate(cluster = "low vs Non DE")
  
  
  dat <- bind_rows(b1, b2, b3, b12, b23, b13) %>% 
    mutate(cluster = factor(cluster, levels = c("low DE", 'high DE', "Non-DE", 
                                                "low vs high DE", "high vs Non DE", "low vs Non DE")))
  
  f1 <- ggplot(data = dat, aes(value, color = cluster)) + 
    geom_density(alpha = 1, aes(fill = cluster), position = "dodge") +
    theme(legend.position = "bottom", 
          axis.title = element_text(face = "bold")) + 
    labs(x = "statitics correlation",
         subtitle = paste0("True data-row correlation is ", rho, ".")) + 
    scale_color_manual(values = c("red", "#00FF00", "#3300FF", 
                                  "#000000", "#33FFFF", "#FF33FF"))
  
  
  return(f1)
}