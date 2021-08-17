//
//  STRATIFIED POWER PRIOR FOR BINARY DATA
//  Power parameter As follows dirichlet
//  Non-centered version. theta = mu + tau * theta_tilda
//
data {
  int<lower = 2> S;

  //existing data
  int<lower = 1> NCH[S];
  int<lower = 0> YSUMCH[S];

  //current data
  int<lower = 1> NE[S];
  int<lower = 1> NCD[S];
  int<lower = 0> YSUME[S];
  int<lower = 0> YSUMCD[S];
  
}

transformed data {
  row_vector<lower = 0, upper = 1>[S] WS1;
  int<lower = 0> sn1;
  
  sn1 = sum(NE);
  for (i in 1:S) {
    WS1[i] = NE[i];
    WS1[i] = WS1[i]/sn1;
  }
}

parameters {
  vector<lower=0, upper=1>[S]  thetaEs;
  vector[S] muCs;
  real <lower=0> tau2s[S];
  real gammaCDs_tilda[S];
  real gammaCHs_tilda[S];
}

transformed parameters {
  vector<lower=0, upper=1>[S]  thetaCs;
  real gammaCDs[S];
  real gammaCHs[S];
  thetaCs = inv_logit(muCs);
  for (i in 1:S){
    gammaCDs[i] = muCs[i] + (sqrt(tau2s[i]))*gammaCDs_tilda[i];
    gammaCHs[i] = muCs[i] + (sqrt(tau2s[i]))*gammaCHs_tilda[i];
  }
}

model {
  //prior
  gammaCDs_tilda ~ normal(0,1);
  gammaCHs_tilda ~ normal(0,1);
  muCs ~ normal(0,10000);
  tau2s ~ inv_gamma(1,0.1);
 
  thetaEs ~ beta(1,1);

  //likelihood
  YSUME ~ binomial(NE, thetaEs);
  YSUMCD ~ binomial_logit(NCD, gammaCDs);
  YSUMCH ~ binomial_logit(NCH, gammaCHs);
}

generated quantities {
  real thetaE;
  real thetaC;
  real delta;
  thetaE = WS1 * thetaEs;
  thetaC = WS1 * thetaCs;
  delta = thetaE - thetaC;
}

