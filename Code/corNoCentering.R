mux <- c(4, 2);
muy <- c(2, 1);


rhovec <- seq(-0.99, 0.99, by = 0.01)
theoretical <- sample <- c()

for ( i in 1:length(rhovec)){
  
  n1 <- 1000;n2 <- 1000; n <- n1 + n2
  rho <- rhovec[i]
  sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  
  library(MASS)
  g1 <- mvrnorm(n=n1, mu = c(mux[1], muy[1]), Sigma = sigma)
  g2 <- mvrnorm(n=n1, mu = c(mux[2], muy[2]), Sigma = sigma)
  
  g <- rbind(g1, g2)
  
  
  d1 <- (rho + n1*n2/n^2*(mux[1]- mux[2])*(muy[1]-muy[2]))
  d2 <- sqrt((1-1/n + n1*n2/n^2*(mux[1]-mux[2])^2)*(1-1/n + n1*n2/n^2*(muy[1]-muy[2])^2))
  theoretical[i] <- d1/d2
  sample[i] <- cor(g)[1, 2]  
  
}

correlation <- data.frame(rhovec, theoretical, sample)
head(correlation)
library(reshape2)
corr1 <- melt(correlation, variable.name = "type", id ="rhovec")
library(ggplot2)
ggplot(data = corr1, aes(x =rhovec, y = value, group = type, color= type)) +
  geom_point() + 
  geom_abline(intercept = 0, slope= 1) + ylim(-1, 1)

