//
//  POWER PRIOR FOR BINARY DATA
//  Power parameter As follows dirichlet
//  Matching + Power Prior. Lin 2019. 
//  Bayes model 1. logit(theta) =beta0 + lambda*Z
//
data {
  int<lower = 1> NECD;
  int<lower = 1> N;
  //data
  real<lower=0> a0[N];
  int Y[N];
  int Z[N];
}

parameters {
  real beta0;
  real lambda;
}

transformed parameters {
  real <lower=0, upper=1> thetaE;
  real <lower=0, upper=1> thetaC;
  real delta;
  thetaE = inv_logit(beta0+lambda);
  thetaC = inv_logit(beta0);
  delta = thetaE-thetaC;
}

model {
  real thetas[N];
  //prior
  beta0 ~ normal(0,10000);
  lambda ~ normal(0,10000);
  
  //likelihood
  for (i in (NECD+1):N) {
    thetas[i]=inv_logit(beta0 + lambda*Z[i]);
    target += bernoulli_lpmf(Y[i]|(thetas[i]))*a0[i];
  }
  for (i in 1:NECD) {
    thetas[i]=inv_logit(beta0 + lambda*Z[i]);
    Y[i] ~ bernoulli(thetas[i]);
  }
}

