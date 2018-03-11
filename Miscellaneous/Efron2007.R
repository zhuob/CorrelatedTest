library(MASS)



##   simple linear regression  correlation of estimate beta
nsim <- 20

x <- cbind( rep(1, 2*nsim), c(rep(0, nsim ), rep(1, nsim)))

beta <- c(10, 2) # gene1  control = 10, trt= 12
mu1 <- t(unique(x%*%beta))

delta <- c(10, 5) # gene2  control = 10 , trt = 15
mu2 <-  t(unique(x%*%delta))

x.inv <- solve(t(x)%*%x)






# generate correlated poisson random variable

cor.pois <- function(n, mu1, mu2, mu3)
{
  y1 <- rpois(n, mu1)
  y2 <- rpois(n, mu2)
  y3 <- rpois(n, mu3)
  
  x1 <- y1  + y2
  x2 <- y1  + y3
  
 # the correlation between x1 and x2 
 # rho <- mu1 /sqrt( (mu1 + mu3)*(mu1 + mu2))
#  print(round(rho, 4))
  return(cbind(x1, x2))
}



t.stat <- matrix(NA, 10000, 3)
z= rho = matrix(NA, 10000, 2)
num <- 20  # how many samples in each group 


for (i in 1:10000)
{
  # in this case, gene 1 is not DE, but gene 2 is DE
  # for treatment 1
  exp1 <- cor.pois(num, 10, 30, 40) # pi1 = 4, p2= 5, rho = 0.2236
  # for treatment 2
  exp2 <- cor.pois(num, 20, 20, 180) #pi1 = 4, p2= 20 rho = 0.2236
  
  t1 <- as.numeric(t.test(exp1[, 1], exp2[, 1])$stat)
  t2 <- as.numeric(t.test(exp1[, 2], exp2[, 2])$stat, alternative="less")
  
  t3 <- as.numeric(t.test(exp1[, 2], exp2[, 2]-150)$stat)
  
  
  z[i,] = colMeans(exp2) - colMeans(exp1)
  
  rho[i, ] <- c(cor(exp1)[1, 2], cor(exp2)[1, 2])
  t.stat[i, ] <- c(t1, t2, t3)
}

cor(t.stat) # correlation of t.stat is the same as cor(zval)
apply(rho, 2, mean)
pval <- pt(t.stat, 2*num-2) 
cor(pval)    #  but not for cor(pval)
zval <- qnorm(pval)
cor(zval)


# 
mu3 <- 3
mu4 <- 4
mu5 <- 6
mu1 <- 0.5*mu4
mu2 <- 0.5*mu4 + mu5
mu6 <- mu4 + 4*mu3

c(mu1, mu2, mu3, mu4, mu5, mu6)
rho1 <- mu1 /sqrt( (mu1 + mu3)*(mu1 + mu2))
rho2 <- mu4 /sqrt( (mu4 + mu5)*(mu4 + mu6))
c(rho1, rho2)



# How close is a t distribution to normal 
df <- 10
va <- df/(df-2)
pi <- seq(-2, 2, length=1000)

f1 <- dnorm(pi, 0, sqrt(va))
f2 <- dt(pi, df)

plot(pi, f1, pch=20, cex=0.3)
lines(pi, f2, pch=3, cex=0.3, col="red")


## how is 77 obtained in page 100 of Efron's paper
N <- 3226
A <- 0.57
x0 <- 2.5
bar.phi <- 1-pnorm(x0)
phi <- dnorm(x0)
cond <- N*bar.phi*(1 + A*x0*phi/(sqrt(2)*bar.phi))
cond

## reproduce figure 2  ------------------------------------------------- #####
##############################################################################

FDP <- FDP_uncondition <- FDP_condition <- c()
source("~/Google Drive/Study/Thesis/Correlation/FDR/locfdr/R/locfns.R")
A <- A.hat1 <- A.hat2 <-c()
sighat1 <- sighat2 <- c()
p.hat <- c()
x <- 2.5
x0 <- 1
COUNT <- 0
for ( i in 1:200){
  
N <- 3000                                                 # number of genes to be simulated
K <- 101                                                  # number of bins

# z1 <- rnorm((0.95)*N, 0, 1)                               # null z values
z2 <- 2.5 + sqrt(1.25)*qnorm(((1:150)-0.5)/150)           # non-null z values
# z <- c(z1, z2)                                            # all z values
bins <- seq(-4.1, 6.0, length=K+1)                        # devide the z values into bins
z.k <- 0.5*(bins[-1] + bins[-K])                          # bin center
z.k[K] <- 5.95
delta <- diff(bins)[1]                                    # bin width

## generating the null counts according to eq(3.15)
alpha <- 0.15                                             # the variance of correlation
A[i] <- rnorm(1, 0, alpha)                                # the dispersion variate (KEY QUANTITY FOR CONDITIONAL FDR)
nu <- N*delta*dnorm(z.k)                                  # the unconditional mean vector for y
W <- N*delta*dnorm(z.k)*(z.k^2-1)/sqrt(2)                 # the wing shape function
u <- nu + A[i]*W                                          # conditional mean given A
u[u<0] <- 0                                               # Truncate u at 0
y1 <- rpois(length(u), u)                                 # eq. 3.15, NULL counts

y2 <- hist(z2, breaks=bins, plot = F)$counts              # nonNULL counts
y <- y1 + y2

## estimating A.hat and sigma0 

# (1) how many bins are in interval [-x0, x0]
y.k <- y[z.k <= 1 & z.k >= -1]
z.k0 <- z.k[z.k <= 1 & z.k >= -1]


## based on my understanding
# model <- glm(y.k ~poly(z.k0, degree=2), family = poisson(link = "log")) 
# beta.2 <- summary(model)$coe[3, 1]                        # eq. 5.16
# sighat1[i] <-  1/sqrt(-2*beta.2)                          # eq. 5.17
# A.hat1[i] <- (sighat1[i]^2 -1)/sqrt(2)
# #sigma.0 <- (1 + sqrt(2)*A[i])^(0.5)

## procedure based on code from locfdr 
library(splines)
f <- glm(y ~ ns(z.k, df= 7), poisson)$fit                  # nonparametric fit
l <- log(f)
imax <- seq(l)[l == max(l)][1]
zmax <- z.k[imax]
x00 <- cbind(z.k0-zmax, (z.k0-zmax)^2) 
y0 <- as.numeric(l[z.k <= 1 & z.k >= -1])
lr <- lm(y0~ x00)
co <- lr$coef
sighat2[i] <- 1/sqrt(-2*co[3])
A.hat2[i] <- (sighat2[i]^2 -1)/sqrt(2)

sigma.0 <- sighat2[i]
## estimating FDRs
P.0.sigma <- 2*pnorm(x0/sigma.0)-1                        # eq. 5.19
Y0 = sum(y[z.k >= -x0 & z.k <= x0])                       # eq. 5.10
p.0.hat <- (Y0/N)/P.0.sigma                               # eq. 5.20
if(p.0.hat >= 1) {COUNT <- COUNT + 1}
  
# X0 <- cbind(1, z.k-zmax, (z.k-zmax)^2)
# p.hat[i] <- sum(exp(as.vector(X0 %*% co)))/sum(f)
#p.0.hat <- p.hat[i]
T.x <- sum(y[z.k > x])                                    # eq. 5.1
FDP_condition[i] <- N*p.0.hat*(1-pnorm(x/sigma.0))/T.x    # eq. 5.21  conditional FDR
FDP_condition[i] <- min(FDP_condition[i], 1)
# print(FDP_condition[i])
p.00 <- (Y0/N)/(2*pnorm(x0)-1 )
FDP_uncondition[i] <- N*p.00*(1-pnorm(x))/T.x             # eq. 5.22 unconditional FDR

FDP[i] <- sum(y1[z.k > x])/T.x                            # eq. 5.3 & 5.5  TRUE FDR

}

# A.hat <- A.hat2
# plot(A, A.hat, cex=0.5)
# abline(lm(A.hat~A))
# abline(h=0, col= "red", lty=2)
# plot(A.hat, FDP_condition, cex=0.5, xlim=c(-0.2, 0.7))
# points(A.hat, FDP_uncondition, cex=0.5, pch=4, col="red")
# srg <- smooth.spline(x=A.hat, y = FDP_condition, df= 3)
# lrg1 <- lm(FDP_condition~A.hat)
# lines(srg, lty=5, lwd=1, col="blue")
# abline(v=0, col="magenta")
# theo.FDR <- 0.95*(1-pnorm(x))/( 0.95*(1-pnorm(x)) + 0.05*0.5)
# abline(h=theo.FDR, col="magenta")


plot(FDP, FDP_condition, cex=0.5, ylim=c(0, 1), xlim=c(0, .6), main=paste("alpha=", alpha))
points(FDP, FDP_uncondition, cex=0.5, pch=4, col="red")
summary(lm(FDP_uncondition~FDP))
lrg <- lm(FDP_condition~FDP)
abline(lm(FDP_uncondition~FDP))
# srg <- smooth.spline(FDP, FDP_condition)
abline(lrg, col="blue", lty=3, lwd=2)



#################################### End ######################################



################# Calculate the quantity of A for HIV data
library(locfdr)
data("hivdata")
z <- hivdata
x0 <- 1
P0  <- 2*pnorm(x0)-1
Q0 <- sqrt(2)*x0*dnorm(x0)
P0.hat <- sum(z <= x0 & z >= -x0)/length(z)
A.hat <- (P0-P0.hat)/Q0

am <- locfdr(z)



##############################################################################

# the z transformation will preserve the correlation structure
z2 <- mvrnorm(100, mu=c(2, 0), sigma)
p <- pnorm(z2)
invz <- qnorm(p)
cor(z2)
cor(invz)



## exploring t.test: 
## the default t.test is two sided
nsim <- 20
x1 <- rnorm(nsim, 0, 1)
x2 <- rnorm(nsim, 1.2, 0.8)
t.test(x1, x2)$p.va

# by hand
s1 <- var(x1)
s2 <- var(x2)
t <- (mean(x1)-mean(x2))/ (sqrt(s1/nsim + s2/nsim))
df <- (s1/nsim + s2/nsim)^2/( 
  (s1/nsim)^2/(nsim-1) + (s2/nsim)^2/(nsim-1))
pt(t, df)*2 



## t to z
# the conclusion here is that, if one of the test is non-null, then
# the correlation is not consistent

nsim <- 100
d <- 20
df <- 20
x <- rt(nsim, df)
z <- qnorm(pt(x + d, df= df))
cor(x, z)

N <- 1000
nsim <- 20
stat.t <- matrix(NA, N, 3)
rho <- matrix(NA, N, 2)

for ( i in 1:N)
{ 
  sig <- 10
  sigma <- sig *matrix(c(1, 0.5, 0.5, 1), 2)
  set.seed(i)
  x1 <- mvrnorm(nsim, mu = c(20, 30), sigma)
  set.seed(1e4 + i+5)
  x2 <- mvrnorm(nsim, mu= c(0, 32), sigma)
  s1 <- apply(x1, 2, var)
  s2 <- apply(x2, 2, var)
  mu1 <- apply(x1, 2, mean)
  mu2 <- apply(x2, 2, mean)
  denom1 <-  (sqrt(s1[1]/nsim + s2[1]/nsim))
  denom2 <- (sqrt(s1[2]/nsim + s2[2]/nsim))
  
 # denom1 <- 1
#    denom2 <- 10
  
  t1 <- (mu1[1]-mu2[1]-18)/denom1
  t2 <- (mu1[2]-mu2[2])/denom2
  t3 <- (mu1[1]-mu2[1])/denom1
  
#   t1 <- t.test(x1[, 1]-18, x2[, 1])$stat
#   t2 <- t.test(x1[, 2], x2[, 2])$stat

    stat.t[i,] <- c(t1, t2, t3)
    rho[i, ] <- c(cor(x1)[1, 2], cor(x2)[1, 2])
}

apply(rho,2, mean)
cor(stat.t)





## plot z.cor, sample.cor and true cor
plot.cor.z.test <- function(mu1, mu2, sig1, sig2, sig3, sig4, sigqua){

true.cor  <- seq(-0.9, 0.9, 0.01) 
z.cor <- samp.cor <- c()

for ( j in 1: length(true.cor))
{
  permu.B <- 100
  rho <- true.cor[j]

# mean 
  #mu1 <- c(20, 30)
  #mu2 <- c(0, 32)

# covariance 
  sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
  sigma2 <- matrix(c(sig3^2, rho*sig3*sig4, rho*sig3*sig4, sig4^2), 2, 2)
  nsim1 <- 20
  nsim2 <- 30

  z.stat <- matrix(NA, permu.B, 2)
  s.cor <- c()
  for ( i in 1:permu.B)
    {
  
      x1 <- mvrnorm(nsim1, mu = mu1, sigma1)
      x2 <- mvrnorm(nsim2, mu = mu2, sigma2)
  
      s1 <- sqrt(sigma1[1,1]/nsim1 + sigma2[1, 1]/nsim2)
      s2 <- sqrt(sigma1[2,2]/nsim1 + sigma2[2, 2]/nsim2)
  
      t1 <- (mean(x1[, 1]) - mean(x2[, 1]))/s1
      t2 <- (mean(x1[, 2]) - mean(x2[, 2]))/s2
  

      
      rm1 <- x1[, 1] - mean(x1[, 1])
      rm2 <- x1[, 2] - mean(x1[, 2])
      
      rm3 <- x2[, 1] - mean(x2[, 1])
      rm4 <- x2[, 2] - mean(x2[, 2])
      s.cor[i] <- cor(c(rm1, rm3), c(rm2, rm4)) # sample corr
      #s.cor[i] <- cor(c(x1[, 1], x2[, 1]), c(x1[, 2], x2[, 2]))
#       p1 <- sum(x1[ ,1]*x1[ ,2])- nsim1*mean(x1[, 1])*mean(x1[, 2])
#       p2 <- sum(x2[ ,1]*x2[ ,2])- nsim2*mean(x2[, 1])*mean(x2[, 2])
#       d1 <- sqrt(sum(x1[, 1]^2) + sum(x2[, 1]^2)- nsim1*mean(x1[, 1])^2 - nsim2*mean(x2[, 1])^2)
#       d2 <- sqrt(sum(x1[, 2]^2) + sum(x2[, 2]^2)- nsim1*mean(x1[, 2])^2 - nsim2*mean(x2[, 2])^2)
#       r <- (p1 + p2)/(d1*d2)
        
      z.stat[i, ] <- c(t1, t2)
      
    }

    z.cor[j] <- cor(z.stat)[1, 2]
    samp.cor[j] <- mean(s.cor)
 

}

plot(true.cor, samp.cor, pch=3, cex=0.4, main=sigqua)
abline(0, 1)
points(true.cor, z.cor, pch=20, cex=0.5, col='red')
legend("topleft", legend=c("sample", "test stat"), pch=c(3, 20), col=c("black", "red"))
}

mu1 <- c(20, 60)
mu2 <- c(5, 30)

plot.cor.z.test(mu1, mu2, sig1=11, sig2= 3,sig3=0.2, sig4=0.2, sigqua="DE")

mu1 <- c(20, 60)
mu2 <- c(200, 120)
plot.cor.z.test(mu1, mu2, sig1=100, sig2= 5,sig3=10, sig4=0.5, sigqua="rho constant")
plot.cor.z.test(mu1, mu2, sig1=10, sig2= 10, sig3=10, sig4=10, sigqua="Null") 
# for constant rho,  sig1*sig4 must equal to sig2*sig3
plot.cor.z.test(mu1, mu2, sig1=0.2, sig2= 0.3,sig3=0.2, sig4=0.3, sigqua="DE")
plot.cor.z.test(mu1, mu2, sig1=50, sig2= 0.3,sig3=0.03, sig4=10, sigqua="DE")



### use pooled correlation? ??
## be cautious 

rho <- 0.3
sig1=50
sig2= 0.1*sig1
mu1 <- mu2 <-  c(0, 0)
sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
rand1 <- mvrnorm(1000, mu1, sigma1)
cor(rand1)

sig4=50
sig3=0.1*sig4

sigma2 <- matrix(c(sig3^2, rho*sig3*sig4, rho*sig3*sig4, sig4^2), 2, 2)
rand2 <- mvrnorm(1000, mu2, sigma2)
cor(rand2)
rand <- rbind(rand1, rand2)
cor(rand)
tt1 <- c(rand1[, 1], rand2[, 1])
tt2 <- c(rand1[, 2], rand2[, 2])
cor(tt1, tt2)





## explore the correlations of two sample t test 
preseve.cor <- function(mu1, mu2, coeff, sig, nsim)  # coeff is the correlation between two genes
{
  
  sigma <- sig*matrix(c(1, coeff, coeff, 1), 2, 2)
  
  #
  t.stat1 <- matrix(NA, nsim, 2) # store the test stats
  # rho1 <- t.stat1 # store the sample corr
  rho1 <- c()
  num1 <- 20
  
  # study correlation 
  for ( i in 1:nsim)
  {
    z1 <- mvrnorm(num1, mu=mu1, sigma)
    #  z1 contains obs from two correlated  genes,  rho=0.6
    z2 <- mvrnorm(num1, mu=mu2, sigma)
    # z2  is just another replicate of the same genes
    #  gene1 Not DE, gene2 DE
    
    # two sample t.test
    t1 <- t.test(z1[, 1], z2[, 1])$stat
    t2 <- t.test(z1[, 2], z2[, 2])$stat
    
    # rho1[i, ] <- c(cor(z1)[1, 2], cor(z2)[1, 2])
    rm1 <- z1[, 1] - mean(z1[, 1])
    rm2 <- z1[, 2] - mean(z1[, 2])
    
    rm3 <- z2[, 1] - mean(z2[, 1])
    rm4 <- z2[, 2] - mean(z2[, 2])
    rho1[i] <- cor(c(rm1, rm3), c(rm2, rm4))
    
    t.stat1[i, ] <- c(t1, t2)
    
  }
  
  t.cor <- cor(t.stat1)[1, 2]# correlation between t test stats
  #p.t <- pt(t.stat1, 2*num1-2) # p value does not preserve the correlation
  #p.cor <- cor(p.t)[1, 2] # correlation between p values
  # mean.coeff <- apply(rho1, 2, mean)
  mean.coeff <- mean(rho1)
  # z <- qnorm(p.t) # z value is similar to test stat and thus sample correlation
  # z.cor <- cor(z)[1, 2] # correlation between z test stats
  
  return(c(t.cor, mean.coeff, coeff))
}



## the correlation cannot be preserved if the null are not true
## even if we consider the test statistics
coeff <- seq(-0.9, 0.9, by = 0.01)
mu1 <- c(20, 30)
mu2 <- c(10, 40)
sig <- 10  # if the variance is too large, the correlation preserves 
corr.matrix <- data.frame(matrix(NA, length(coeff), 3))
colnames(corr.matrix) <- c("t","rho", "rho.true")

for (i in 1:length(coeff))
{
  corr.matrix[i, ] <- preseve.cor(mu1, mu2, coeff[i], sig, nsim = 100)
  
}

plot(corr.matrix$rho.true, corr.matrix$rho, pch=20, cex=0.3, 
     xlab="rho.true", ylab="rho",
     main = paste("mu1 = ", mu1, "mu2 = ", mu2, "sigma =", sig))
abline(0, 1)
points(corr.matrix$rho.true, corr.matrix$t, cex= 0.5, pch =3, col="magenta")


s1 <- matrix(NA, 10000, 2)

for ( i in 1: 10000)
{
  
  rho <- 0.6
  sig1=50
  sig2= 0.8*sig1
  mu1 <- mu2 <-  c(0, 0)
  sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
  rand1 <- mvrnorm(1000, mu1, sigma1)
  v <- apply(rand1, 2, var)
  s0 <- v[1]/sig1^2
  s2 <- v[2]/sig2^2
  s1[i, ] <- c(s0, s2)
}
cor(s1)




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





########### Verify that the correlation has mean 0 --------------------------
########### Explore the variance of correlation if DE is introduced ----------

size <- 50                                     # number of samples to be considered
n_gene <- 500                                  # number of genes to be considered
de_prop <- 0.0                                 # proportion of DE
de <- n_gene*de_prop                           # number of DE genes
delta <- rnorm(n_gene, 0, 4)                     # magnitude of DE
z <- rbinom(n_gene, 1, de_prop)                # indicator, DE =1, NON_DE=0 
delta <- delta*z
#delta <- rep(0, n_gene)

rho <- 0.3                                     # correlation matrix
sigma0 <- 1                                    # variance part 
cor.struct <- matrix(rho, n_gene, n_gene);    
diag(cor.struct) <- 1                          # correlation matrix
sigma <- sigma0*cor.struct                     # covariance matrix for both groups

mu_initial <- 0
mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
mu2 <- delta  + mu_initial                     # mean expression value for treatment

library(MASS)
microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix

n <- size/2
#set.seed(50)
microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
#set.seed(51)
microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))


##  Convolution to get the true correlation matrix 
sigma <- sigma0*cor(t(microarray))
microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))

sigma <- sigma0*cor(t(microarray))
microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))

sigma <- sigma0*cor(t(microarray))
microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))


microarray_scale <- scale(microarray)         # standardize the microarray LIKE NORMALIZATION

####  Nullify any treatment effect
source("SimulateLabData.R")
trt <- rep(c(0, 1), each= size/2)                              # group indicator
group_mean <- as.matrix(group.mean(microarray_scale, trt))     # calculate correlation matrix
resid_mat <- microarray_scale - group_mean                     # the trt effects are removed from matrix

samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
c(mean(samp_rho), sd(samp_rho))
hist(samp_rho)





####################################################################################
#    create an example where correlation and nonnull have the same hist of z values
#    see notes on Week 0 of 2015 Fall
####################################################################################

# Case 1: z[1]-z[900] ~ N(0, 1) as null parts, z[901]--z[1000] ~ N(d, 1), d ~ N(0, 3)
# Note, if z|d ~ N(d, 1) and d ~ N(0, n) then z ~ N(0, n + 1)

# Case 2:  z[1]- z[1000] ~ N(0, 1) but they are somehow correlated.



n <- 1e3                                                 # number of z values
m <- 10
set.seed(100)
z <- rnorm(n*m)
dim(z) <- c(n, m)
z.cov <- cor(t(z))

mu <- rep(0, n)
library(MASS)
set.seed(329)
z.null <- mvrnorm(1, mu, z.cov)
hist(z.null, breaks = seq(-4, 4, by= 0.2), freq = F, main="correlated z, n=1000", ylim=c(0, 0.4))
x <- dnorm(seq(-4, 4, by =0.01))
lines(seq(-4, 4, by = 0.01), x, col="red")

z2 <- rnorm(n)
hist(z2, breaks = seq(-4, 4, by= 0.2), freq = F, main="z ~ N(0, 1), n=1000")
x <- dnorm(seq(-4, 4, by =0.01))
lines(seq(-4, 4, by = 0.01), x, col="red")

seed <- 32
p0 <- 0.90
z <- rep(NA, n)
n.null <- p0*n
set.seed(101)
z[1:(p0*n)] <- rnorm(n.null)
#set.seed(seed)
delta <- rnorm(n-n.null, 0, 2)
z[-(1:(p0*n))] <- rnorm(n-n.null, delta, 1 )

z <- z[z>-4 & z < 4]
hist(z, breaks = seq(-4, 4, by= 0.2), freq = F, main="mixed normal, n=1000")
x <- dnorm(seq(-4, 4, by =0.01))
lines(seq(-4, 4, by = 0.01), x, col="red")

# hist(z.null, freq=F, breaks=30)
# x <- seq(-4, 4, by = 0.01)
# lines(x, dnorm(x), col= "red")
# # c(mean(z.correlated), var(z.correlated))
# p = pmin(pnorm(z.null), pnorm(z.null, lower.tail=FALSE))*2;
# hist(p);


p0 <- 0.90
#mu.z <- mean(z.null)/(1-p0)
#sigma.z <- (var(z.null)- p0^2)/(1-p0)^2
#c(mu.z, sigma.z)



hseed <- function(seed)
{
  
z <- rep(NA, n)
  n.null <- p0*n
  set.seed(101)
  z[1:(p0*n)] <- rnorm(n.null)
  set.seed(seed)
  delta <- rnorm(n-n.null, 0, 2)
  z[-(1:(p0*n))] <- rnorm(n-n.null, delta, 1 )


### even if z values are from the same distribution, the histogram does not
### necessarily look the same. 
# z1 <- rnorm(1000)
# z2 <- rnorm(1000)
# z.val <- data.frame(test =c(z1, z2), group = rep(c("True Alternative", "Correlated Z"), each = n))

z.val <- data.frame(test =c(z, z.null), group = rep(c("True Alternative", "Correlated Z"), each = n))
  
  
#now make your lovely plot
library(ggplot2)
ggplot(z.val, aes(test, fill = group)) + geom_density(alpha = 0.2,aes(y = ..density..)) + 
  stat_function(fun = dnorm, colour = "red")  +                  ## add a normal density
  xlim(-6, 6)

}


hseed(215)
hseed(222)
hseed(233)








##############################################################################
## column standardization makes sum of the correlation 0 
##############################################################################

r.mat <- matrix(.3, 500, 500)
diag(r.mat) <- 1
library(MASS)
n <- 15
x <- matrix(NA, 500, n)

for (k in 1: n)
{
  #mu <- rep(k, 500)
  if (k < n/2)  {mu <- rnorm(500)}
  else { mu <- rnorm(500, 1, 1) }
  xk <- mvrnorm(n=1, mu, Sigma = 2*r.mat)                 # multivariate random number
  x[, k] <- (xk)                                     # column standardization
}

x.standard <- scale(x)

# rho <- cor(t(x))
rho <- cor(t(x.standard))
hist(rho, breaks=30, freq=F, main=paste("mean =", round(mean(rho), 3)))



############################################################################################
##########################  after standardization, there's no treatment effect ########
############################################################################################
p1 <- p2 <- c()
for ( i in 1: nrow(x))
{
  p1[i] <- t.test(x[i, 1:floor(n/2)], x[i, -(1:floor(n/2))])$p.val
  p2[i] <- t.test(x.standard[i, 1:floor(n/2)], x.standard[i, -(1:floor(n/2))])$p.val
  
}

hist(p1, breaks=30)
hist(p2, breaks=30)

















######################  TAIL Counts #################

library(MASS);

m = 1000;
n = 10;

histz = function(seed) {
  set.seed(seed);
  z = rnorm(m*n);
  dim(z) = c(m,n);
  
  v = cor(t(z));
  
  mean(v);
  
  det(v);
  ## hist(v);
  quantile(v);
  
  mu = rep(0,m);
  
  z = mvrnorm(10, mu, v);
  
  hist(z, prob=TRUE, breaks=30);
  x = seq(-4,4,length=100);
  lines(x, dnorm(x));
  
  p = pmin(pnorm(z), pnorm(z, lower.tail=FALSE))*2;
  hist(p);
  
}

par(mfcol=c(4,4));
## for (i in 1001:1016){
for (i in 1001:1008){
  histz(i);
}

for (i in 1009:1016){
  histz(i);
}

par(mfrow=c(1, 2))
histz(999);

histz(1001);

## y = apply(z, 2, function(x) scale(x));
## v = cor(t(y));

## simulation of P(z1 > z0, z2 > z0) when (z1, z2) are bivariate normal.
pz0 <- function(rho, z0)
{
  mu <- c(0, 0)
  sigma <- c(1, rho, rho, 1)
  dim(sigma) <- c(2, 2)
  
  x <- mvrnorm(n=1000, mu, sigma)
  index <- x[, 1] > z0 & x[, 2] > z0
  p0 <- mean(index)                  
  return(c(rho, p0))
}

rho <- seq(-0.99, 0.99, by=0.01)
p0 <- c()
par(mfrow=c(1, 1))
z0 <- seq(0.5, 3, by=0.1)
z0 <- 2
for( i in 1:length(rho))
{p0[i] <- pz0(rho[i], z0)[2]}
plot(rho, p0, cex=0.5, pch=20)




##### for locfdr package, all that matters is just the z-vlaues, .....
##### through which Efron calculated A and sigma and conditional FDR.....
##### Basically, that means if I can generate a vector of z-vlaues with ....
##### the same histogram, this procedure cannot tell whether true DE or correlation contributes....

n <- 1e3                                                 # number of z values
m <- 10
z <- rnorm(n*m)
dim(z) <- c(n, m)
z.cov <- cor(t(z))

mu <- rep(0, n)
library(MASS)
#set.seed(299)
z.null <- mvrnorm(1, mu, z.cov)
par(mfrow=c(1, 2))
h <- hist(z.null, freq=F, breaks=30, ylim=c(0, 0.5))
x <- seq(-4, 4, by = 0.01)
lines(x, dnorm(x), col= "red")
p = pmin(pnorm(z.null), pnorm(z.null, lower.tail=FALSE))*2;
hist(p);
c(mean(z.null), sd(z.null))

library(locfdr)

fdr <- locfdr(z.null, bre=30)
probs <- h$counts/n   # probability in each bin
count <- r








