//
//  IPTW + Commensurate prior (BINARY DATA)\
//  Bayes model 1. logit(theta) =b0 + b1*Z
//  Non-centered version. theta = mu + tau * theta_tilda
//  theta_CD ~ N(theta_CH, sigma2)
//
data {
  // Prior parameters. tau2 ~ inv_gamma(prior_alpha=1, prior_beta=0.1);
  real<lower = 0> prior_alpha;
  real<lower = 0> prior_beta;  
  int<lower = 1> NE;
  int<lower = 1> N;
  int<lower = 1> NCD;
  real<lower=0> a0[N];
  //data
  int Y[N];
}

parameters {
  real <lower=0, upper=1> thetaE;
  real gammaCH;
  real gammaCD;
  real <lower=0> tau2;
}

transformed parameters {
  real delta;
  real <lower=0, upper=1> thetaC;
  real <lower=0, upper=1> thetaCH;
  
  thetaC = inv_logit(gammaCD);
  delta = thetaE-inv_logit(gammaCD);
  thetaCH = inv_logit(gammaCH);
}

model {
  //prior
  thetaE ~ beta(1,1);
  gammaCH ~ normal(0,10000);
  tau2 ~ inv_gamma(prior_alpha,prior_beta);
  gammaCD ~ normal(gammaCH,tau2);
  
  //likelihood
  for (i in 1:NCD) {
    Y[i] ~ bernoulli_logit(gammaCD);
  }
  for (i in (1+NCD):(NE+NCD)) {
    Y[i] ~ bernoulli(thetaE);
  }
  for (i in (1+NE+NCD):(N)) {
    target += bernoulli_lpmf(Y[i]|(thetaCH))*a0[i];
  }
}

