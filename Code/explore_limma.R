library(dplyr)
library(limma)
library(ggplot2)
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/product.moments.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/plot.rhoT.vs.rho.R")
source('~/Desktop/correlated/limma-functions.R')


## limma simulation 
# 1000 genes, sample size 10 x 2, 10000 data sets
# in each data set, 5% low DE, 5% high DE, 90% non DE 
# the data row have pairwise correlation of 0.1 or 0.5 (choose wisely)
#    a. for each data set, simulate the genes, run limma, get the 1000 values
#    b. repeat a 10000 times
#    c. compute the 1000x1000 data correlation matrix

nsim <- 1e4; # nubmer of simulations
seed <- 123; # set seed for reproducibility
n1 <- n2 <- 10; # sample size for each group
de_low <- -3  # DE effect for low DE genes
de_high <- 2  # DE effect for high DE genes
n_gene <- 1e3  # number of genes simulated


rho <- 0.5
delta <- c(de_low, de_high, 0)
params <- expand.grid(deltaX = delta, deltaY = delta)
## theoretical correlation under two-sample t-test
if(n1 > 100){
  rho_t <- purrr::map2_dbl(.x = params$deltaX, .y = params$deltaY, .f = compute.rhoTInf, 
                           rho = rho) %>% unique()
} else{
  rho_t <- purrr::map2_dbl(.x = params$deltaX, .y = params$deltaY, .f = compute.rhoT, 
                           v = n1 + n2 -2, rho = rho, n1 = n1, n2 = n2) %>% unique()
}

sim1 <- simulate_t(nsim = nsim, seed = seed, rho = rho, n1 = n1, n2 = n2,
                   de_low = de_low, de_high = de_high, n_gene = n_gene)
# saveRDS(sim1, "sim1.RDS")
# sim1 <- readRDS("sim1.RDS")
rs11 <- plot_t_cor(sim1$t_limma, rho = rho)
rs12 <- plot_t_cor(sim1$tt, rho = rho)
postscript("limma1.eps", width = 6, height = 6)
rs11 + geom_vline(xintercept = rho, linetype = "dotted", color = "black") + 
  labs(title = "Density plot of moderated t-statistics correlation",
       caption = "Vertical line indicates theoretical data-row correlation.")
dev.off()
rs12 + geom_vline(xintercept = rho_t, linetype = "dotted") + 
  labs(title = "Density plot of regular t-statistics correlation",
       caption = "Vertical lines indicate theoretical value of test statistic correlation.")



rho <- 0.1

if(n1 > 100){
  rho_t <- purrr::map2_dbl(.x = params$deltaX, .y = params$deltaY, .f = compute.rhoTInf, 
                           rho = rho) %>% unique()
} else{
  rho_t <- purrr::map2_dbl(.x = params$deltaX, .y = params$deltaY, .f = compute.rhoT, 
                           v = n1 + n2 -2, rho = rho, n1 = n1, n2 = n2) %>% unique()
}

sim2 <- simulate_t(nsim = nsim, seed = seed, rho = rho, n1 = n1, n2 = n2,
                   de_low = de_low, de_high = de_high, n_gene = n_gene)
# saveRDS(sim2, "sim2.RDS")
# sim2 <- readRDS("sim2.RDS")
rs21 <- plot_t_cor(sim2$t_limma, rho = rho)
rs22 <- plot_t_cor(sim2$tt, rho = rho)
postscript("limma2.eps", width = 6, height = 6)
rs21 +  geom_vline(xintercept = rho, linetype = "dotted", color = "black") + 
  labs(title = "Density plot of moderated t-statistics correlation",
       caption = "Vertical line indicates theoretical data-row correlation.") 
dev.off()
rs22 + geom_vline(xintercept = rho_t, linetype = "dotted") + 
  labs(title = "Density plot of regular t-statistics correlation",
       caption = "Vertical lines indicate theoretical value of test statistic correlation.")
