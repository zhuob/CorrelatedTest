## Mon Mar  5 11:14:32 PST 2018

## 1. Implement Joarder's formula for computing product-moments as an
## infinite sum.
##
## 2. Compute Cv(rho) on a set of (v, rho) values and save the results
## to a file.

## Joader (2009) derived a formula for computing 
##
##   E(U1^a U2^b)
##
## as an infinite sum, where (U1, U2) = (v S_X^2/sigma_X^2, vS_y^2/sigma_Y).
##
## Note that we need to be careful to avoid overflow when computing the sum: e.g., converting some terms to log.
##
## The cases when rho = 0 or 1 need special treatment.


## Compute  E(U1^a U2^b)
euab = function(a, b, m, rho, K=10000) {

  rho = abs(rho);
  
  if (rho==0) {
    return(m*2^(a+b)/gamma(m/2)^2*gamma(m/2+a)*gamma(m/2+b));
  }

  if (rho==1) {
    return(m/(m-2));
  }
 
  g= function(k) {
    ##   print(c(lgamma((k+m)/2+a), lgamma((k+m)/2),lgamma((k+m)/2+b),lgamma(m/2),lgamma((k+1)/2),lgamma(k+1)));
    lgamma((k+m)/2+a)- lgamma((k+m)/2)+lgamma((k+m)/2+b)-lgamma(m/2)+lgamma((k+1)/2)-lgamma(k+1);
  }

  
  k = (K:0)*2;
  gk = g(k) + k*log(2*rho)+(a+b+m/2)*log(1-rho^2);
  
  m*2^(a+b)/sqrt(pi)*sum(exp(gk));

}

test.euab = function() {
  rm(list=ls());
  source('product.moments.R');

  euab(-1/2, -1/2, 100, 0);
  euab(-1/2, -1/2, 100, 1e-8);

  a = b = -1/2;
  m = 58;
  rho = 0.05;
  K = 100;
}

Cvrho = function(v, rho) {
  euab(-1/2, -1/2, v, rho);
}

main = function() {
  rm(list=ls());
  source('product.moments.R');
    
  v = c(4, seq(8, 98, 10));
  nv = length(v);
  
  rho = (0:20)*0.05;
  nrho = length(rho);

  cvrho = matrix(0, nv, nrho);
  rownames(cvrho) = v;
  colnames(cvrho) = rho;

  for (i in 1:nv) {
    for (j in 1:nrho) {
      cvrho[i, j]= Cvrho(v[i], rho[j]);
    }
  }

  write.csv(file="cvrho.csv", cvrho);

  matplot(rho, t(cvrho), type="l");

}


