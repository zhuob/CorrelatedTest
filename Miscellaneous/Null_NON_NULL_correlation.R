

#########
### my point: null/nonnull affects the location of mean(z values)
### while correlation affects the variance of z values

##  if the correlation increases, the variance of z decreases.
library(MASS)
Cov.1 <- c()

rho <- seq(0.1, 0.99, by = 0.01)
for ( k in 1:length(rho))
{
  ndim <- 100
  mu <- rep(0, ndim)
  V <- matrix(rho[k], ndim, ndim)
  diag(V) <- 1
  sigma <- V * 1
  rand.effect <- mvrnorm(1, mu, V)
  Cov.1[k] <- var(rand.effect)
  
}

plot(rho, Cov.1, pch=3, cex=0.5)

abline(lm(Cov.1 ~rho))

### under the null, correlation affects enrichment analysis
rho <- 0.8
ndim <- 200

mu <- rep(0, ndim)
V <- matrix(rho, ndim, ndim)
diag(V) <- 1
sigma <- V*1

correlated.z <- mvrnorm(1, mu, V)
uncor.z <- rnorm(ndim)

beta0 <- 0
pi.uncor <- exp(beta0 + uncor.z)/( 1 + exp(beta0 + uncor.z) )
pi.cor <- exp(beta0 + correlated.z) / ( 1 + exp(beta0 + correlated.z))

GO.uncor <- rbinom(ndim, size=1, pi.uncor)
sum(GO.uncor)
GO.cor <- rbinom(ndim, size=1, pi.cor)
sum(GO.cor)


x <- rnorm(100,2,1)
y <- 1 + 2*x + rand.effect
y2 <- 1 + 2*x + rnorm(100)

summary(lm(y2~ x))
