## Sun Mar 11 20:47:20 PDT 2018
## Add simualted rhoT

## Mon Mar  5 11:56:31 PST 2018
#source("product.moments.R");
#source("simulate.t.test.R");


Cvrho.vec = function(v, rho.vec) {
  unlist(lapply(rho.vec, function(rho){Cvrho(v, rho)}));
}

test.Cvrho.vec = function() {
  source('plot.rhoT.vs.rho.R');
  v = 4;
  rho = seq(-1, 1, 0.05);
  Cvrho.vec(v, rho);
}


compute.rhoT = function(v, rho, deltaX, deltaY, n1, n2) {
  beta = 1/(v * (1/n1 + 1/n2));
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

 


