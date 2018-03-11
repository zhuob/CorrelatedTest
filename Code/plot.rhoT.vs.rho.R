## Mon Mar  5 11:56:31 PST 2018
source("product.moments.R");


Cvrho.vec = function(v, rho.vec) {
  unlist(lapply(rho.vec, function(rho){Cvrho(v, rho)}));
}

test.Cvrho.vec = function() {
  source('plot.rhoT.vs.rho.R');
  v = 4;
  rho = seq(-1, 1, 0.05);
  Cvrho.vec(v, rho);

}


A = function(v) {
  v/(v-2);
}

## E(S_X^-1)=B(v)\sigma_X^-1
B = function(v) {
  sqrt(v/2)*gamma(v/2-1/2)/gamma(v/2);
}

compute.rhoT = function(v, rho, deltaX, deltaY, beta=1/4) {
  Cv = Cvrho.vec(v, rho);
  Bv = B(v);
  Av = A(v);
  r =  v*(Cv - Bv^2);

  rhoT1 = rho * Cv  + beta * deltaX * deltaY * r;
  rhoT2 = Av + beta * v * (Av - Bv^2) * deltaX^2;
  rhoT3 = Av + beta * v * (Av - Bv^2) * deltaY^2;
  rhoT = rhoT1/sqrt(rhoT2 * rhoT3);

}

compute.rhoTInf = function(rho, deltaX, deltaY, beta=1/4) {
  rhoTInf1 = rho  + beta * deltaX * deltaY * 0.5 * rho^2;
  rhoTInf2 = 1 + beta * 0.5 * deltaX^2;
  rhoTInf3 = 1 + beta * 0.5 * deltaY^2;
  rhoTInf = rhoTInf1/sqrt(rhoTInf2 * rhoTInf3);

}

 
plot.rhoT.vs.rho = function(deltaX, deltaY, v) {
  nv = length(v);
  
  rho = seq(-1, 1, length=81);
  nrho = length(rho);

  rhoT = matrix(0, nv+1, nrho);
  rownames(rhoT) = c(v,Inf);
  colnames(rhoT) = rho;
  
  rhoTInf = compute.rhoTInf(rho, deltaX, deltaY);

  rhoT[nv+1, ] = rhoTInf;

  plot(rho, rhoTInf, type="l", lwd=2, col="magenta",
       main=sprintf("(%d, %d)", deltaX, deltaY),
       ylab="rho_T", ylim=c(-1,1));

  for (i in 1:nv) {
    rhoT[i,] = compute.rhoT(v[i], rho, deltaX, deltaY);
    lines(rho, rhoT[i,]);
  }

  list(rho, rhoT);
  abline(0,1);
}

