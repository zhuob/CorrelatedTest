# OS<-.Platform$OS.type

generate.pos <- function(mu, rho, nvectors){  # generate correlated poisson random varaible with
                                           # mu as mean vector of length 2, rho the correlation. 
  
# This is the bivariate normal expression in equation (7) of Madsen and Birkes. The argument c
# corresponds to the parameter delta in (7). If you change the target marginal distribution,
# replace 'ppois' with the new target pmf. If your target marginals require more parameters
# they must be added to the argument list.
g<-function(y,mu1,mu2,c) {
  Sigma<-matrix(c(1,c,c,1),2,2)  
  g<-(pmvnorm(upper=c(qnorm(ppois(y[[1]],mu1)),qnorm(ppois(y[[2]],mu2))),sigma=Sigma)[1])
}

   library(mvtnorm)
  filename<-'uniquecombSAWindows.txt' # File to hold results of numerical solutions.
# f is the function used to numerically solve equation (7) of Madsen and Birkes. The argument c
# is "delta" in equation (7). The argument rho_sb is the target Pearson correlation on the left-hand
# side of equation (7). If you change the target marginals, change 'qpois' to the new marginal cdf
# and 'ppois' to the new marginal pmf. If your target marginals require more parameters
# they must be added to the argument list.
 f<-function(c,mu1,mu2,rho_sb) {
    t<-4 # This is just a "tuning" parameter to get big enough finite upper limits for sums in (7).
    K1<-t*ceiling(qpois(sqrt(1-.005),mu1)) # K1 and K2 are upper limits for sums to infinity in (7).
    K2<-t*ceiling(qpois(sqrt(1-.005),mu2))-1 
    s<-0   # s will hold the value of the double sum on the right-hand-side of (7).
    for (i in seq(0,K2-1,t)) {
      y<-expand.grid(x=seq(0,K1),y=seq(i,i+t-1,1))
      ylist<-split(y,as.numeric(row.names(y)))
      # The function g in this next line is the bivariate normal cdf expression in (7) and is
      # defined in the code above.
      s<-s+sum(1-ppois(y[,'x'],mu1)-ppois(y[,'y'],mu2)+unlist(lapply(ylist,g,mu1=mu1,mu2=mu2,c=c)))
    }
    (s-mu1*mu2)/sqrt(mu1*mu2)
    return((s-mu1*mu2)/sqrt(mu1*mu2)-rho_sb) # From equation (7) in Madsen and Birkes.
  }

# Set parameters for the simulation.
# rho<-.5
# n<-2
# t<-1:n
# muY<-c(2, 3)
# varY<-muY
n <- 2
rho <- rho
muY <- mu  # mean vector
VarY <- muY
 
# CorrY target correlation matrix
# CorrY <- diag(n)
# CorrY <- rho^abs(row(CorrY)-col(CorrY))

CorrY <- matrix(c(1, rho, rho, 1), 2, 2) # correlation matrix

# This section figures out all the unique combinations of pairs of means and correlations
# so we don't have to solve the equation (7) of Madsen & Birkes any more times than necessary.
# If your target marginal has more than just one parameter (like the Poisson), you'll need to include
# those extra parameters here.
mu1<-kronecker(matrix(1,1,n),muY)[upper.tri(kronecker(matrix(1,1,n),muY))]
mu2<-t(kronecker(matrix(1,1,n),muY))[upper.tri(t(kronecker(matrix(1,1,n),muY)))]
c<-CorrY[upper.tri(CorrY)]
c.star<-as.numeric(NA_character_)
mvc<-cbind(mu1,mu2,c,c.star)
unique.comb<-data.frame(unique(mvc))
n.unique<-dim(unique.comb)[1]

# This section goes through each unique combination of marginal parameters and target correlation,
# solves equation (7) and records 'c.star', the 'delta' that solves equation (7). Each line is written
# to the file 'filename'. If you're not getting solutions from uniroot, try just running uniroot for a
# particular i, instead of enclosing it in try(...,silent=TRUE). Also remember to add extra marginal
# parameters to the call to uniroot.
ptm<-proc.time() 
for (i in 1:n.unique) {
  if (unique.comb[i,]$c==0) {
    unique.comb[i,]$c.star<-0
  } else { 
    try(
      unique.comb[i,]$c.star<-uniroot(f,c(max(-1,unique.comb[i,]$c-.1),min(1,unique.comb[i,]$c+.1)),tol=.001,mu1=unique.comb[i,]$mu1,mu2=unique.comb[i,]$mu2,rho_sb=unique.comb[i,]$c)$root
      ,silent=TRUE) 
  }
  write(c(unique.comb[i,]$mu1,unique.comb[i,]$mu2,unique.comb[i,]$c,unique.comb[i,]$c.star),filename, ncolumns = 4, append = TRUE, sep = "\t")  
}
et<-proc.time() - ptm
n.unique

# This section takes the solutions of equation (7) and populates the matrix mvc, 
# which has a row for each pair of points. If your target marginals have more parameters than the
# Poisson, note that you'll need to account for them here.
# These two lines are in case we need to reread the unique.comb data frame:
# unique.comb<-read.table(filename,header=FALSE)
# names(unique.comb)<-names(as.data.frame(mvc))
rlabels<-kronecker(matrix(1,1,n),1:n)[upper.tri(kronecker(matrix(1,1,n),1:n))]
clabels<-t(kronecker(matrix(1,1,n),1:n))[upper.tri(t(kronecker(matrix(1,1,n),1:n)))]
for (i in 1:n.unique) {
  # If target marginals have more than one parameter each, change "1:3" in the line below accordingly.
  # It's "1:3" for Poisson marginals because each marginal has one parameter and they share a correlation parameter.
  mvc[which(abs(rowSums(kronecker(as.matrix(unique.comb[i,1:3]),matrix(1,nrow(mvc),1))-mvc[,1:3]))<=.Machine$double.eps^.25),'c.star']<-unique.comb[i,'c.star']
}
 
# This part does the simulation.
# ptm<-proc.time() 
# nvectors<-1000 # Number of datasets to simulate. 

CorrZ<-matrix(0,n,n) # CorrZ will be the correlation matrix of the standard normal vector.
CorrZ[cbind(rlabels,clabels)]<-mvc[,'c.star']  # Turn the c.star column of mvc into the matrix CorrZ. At first, it's upper triangular
CorrZ<-CorrZ+t(CorrZ)          # Get the lower triangle too.
diag(CorrZ)<-1
R<-chol(CorrZ) # To simulate correlated standard normals.
Y<-array(0,c(n,nvectors)) # Initialize Y to all zeros. The array Y will contain the simulated data.
U<-Y
Kp<-ceiling(max(10*sqrt(muY))) # Kp is the biggest value of Y that will be simulated.
# In the next line, change 'ppois' to the target pmf.
probs<-matrix(unlist(lapply(as.list(0:Kp),ppois,muY)),length(muY),Kp+1) # probs[i,j] is P(Y_i=j-1)
for (v in 1:nvectors) {
	U[,v]<-pnorm(t(R)%*%rnorm(n))  #The U's are uniform(0,1).
	for (j in 1:Kp) {
    Y[which(probs[,j]<U[,v] & U[,v]<=probs[,j+1]),v]<-j
	}
}
et1<-proc.time() - ptm

return(Y)

}
# 
# # Write the simulated data to a file.
# write.table(Y,"YSim.txt",row.names=FALSE,col.names=FALSE)
# 
# #Y<-read.table("YSim.txt",header=FALSE)
# #
# # This next section produces plots to compare sample correlations, means, and variances vs. targets.
# win.graph() 
# ObsCorr<-cor(t(Y))
# xmax<-1.1*max(ObsCorr[lower.tri(ObsCorr)],CorrY[lower.tri(CorrY)])
# xmin<-min(ObsCorr[lower.tri(ObsCorr)],CorrY[lower.tri(CorrY)])-0.1*min(ObsCorr[lower.tri(ObsCorr)],CorrY[lower.tri(CorrY)])
# plot(CorrY[lower.tri(ObsCorr)],ObsCorr[lower.tri(ObsCorr)],main='Figure 1(c)',xlim=c(xmin,xmax),ylim=c(xmin,xmax),xlab="Target Pearson Correlation",ylab="Sample Pearson Correlation",cex.axis=1.5,cex.lab=1.5)
# lines(c(xmin,xmax),c(xmin,xmax),lty=2)
# 
# win.graph()
# ObsMean<-rowMeans(Y)
# xmax<-1.1*max(ObsMean,muY)
# xmin<-xmax/1.1-xmax
# plot(muY,ObsMean,main='Figure 1(a)',xlab="Target Means",ylab="Sample Means",cex.axis=1.5,cex.lab=1.5)
# lines(c(xmin,xmax),c(xmin,xmax),lty=2)
# 
# win.graph()
# ObsVar<-diag(cov(t(Y)))
# xmax<-1.1*max(ObsVar,varY)
# xmin<-xmax/1.1-xmax
# plot(varY,ObsVar,main='Figure 1(b)',xlab="Target Variances",ylab="Sample Variances",cex.axis=1.5,cex.lab=1.5)
# lines(c(xmin,xmax),c(xmin,xmax),lty=2)
