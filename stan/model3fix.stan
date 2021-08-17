//
//  STRATIFIED POWER PRIOR FOR BINARY DATA
//  Fixed power prior
//
data {
  int<lower = 2> S;

  //existing data
  int<lower = 1> NCH[S];
  real<lower = 0, upper = 1> YBARCH[S];

  //current data
  int<lower = 1> NE[S];
  int<lower = 1> NCD[S];
  int<lower = 0> YSUME[S];
  int<lower = 0> YSUMCD[S];

  //prior of vs
  vector<lower=0>[S] RS;

  //target borrowing
  real<lower = 0> A;
}

transformed data {
  row_vector<lower = 0, upper = 1>[S] WS1;
  int<lower = 0> sn1;
  real srs;
  simplex[S] vs;
  real<lower = 0, upper = 1> as[S];
  real<lower = 0> alphas[S];
  real<lower = 0> betas[S];

  sn1 = sum(NE);
  srs = sum(RS);
  for (i in 1:S) {
    WS1[i] = NE[i];
    WS1[i] = WS1[i]/sn1;
    vs[i] = RS[i]/srs;
    as[i]  = 1 < A*vs[i]/NCH[i] ? 1:A*vs[i]/NCH[i];//get min(A*vs[i]/NCH[i], 1)
    alphas[i] = as[i] * NCH[i] * YBARCH[i]  + 1;
    betas[i]  = as[i] * NCH[i] * (1-YBARCH[i]) + 1;
  }
}

parameters {
  vector<lower=0, upper=1>[S]  thetaEs;
  vector<lower=0, upper=1>[S]  thetaCs;
}

model {
  //prior
  target += beta_lpdf(thetaCs | alphas, betas);
 for (i in 1:S) {
     thetaEs[i] ~ beta(1,1);
   }

  //likelihood
  YSUME ~ binomial(NE, thetaEs);
  YSUMCD ~ binomial(NCD, thetaCs);
}

generated quantities {
  real thetaE;
  real thetaC;
  real delta;
  thetaE = WS1 * thetaEs;
  thetaC = WS1 * thetaCs;
  delta = thetaE - thetaC;
}

