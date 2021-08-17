//
//  POWER PRIOR FOR BINARY DATA
//  Power parameter As follows dirichlet
//  Matching + Power Prior. Lin 2019. 
//  Bayes model 2. theta ~ Beta(k*mu, k*(1-mu)) kE=2,muE=0.5
//  Non-centered version. kappC = NC*kappaC_tilda
//
data {
  int<lower = 1> NECD;
  int<lower = 1> NC;
  int<lower = 1> N;
  //data
  real<lower=0, upper=1> a0[N];
  int Y[N];
  int Z[N];
}

parameters {
  real <lower=0> kappaC_tilda;
  real muC;
  real <lower=0, upper=1> thetaE;
  real <lower=0, upper=1> thetaC;
}

transformed parameters {
  real delta;
  real kappaC;
  delta = thetaE-thetaC;
  kappaC = NC*kappaC_tilda;
}

model {
  //prior
  thetaE ~ beta(1,1);
  thetaC ~ beta(kappaC*muC,kappaC*(1-muC));
  kappaC_tilda ~ uniform(0.01,1);
  muC ~ beta(1,1);
  
  //likelihood
  for (i in (NECD+1):N) {
    target += bernoulli_lpmf(Y[i]|(thetaC))*a0[i];
  }
  for (i in 1:NECD) {
    Y[i] ~ bernoulli(thetaC*(1-Z[i])+thetaE*Z[i]);
  }
}

