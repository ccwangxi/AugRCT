########################################
# Title: Compile STAN models (Hybrid methods)
# Author: Xi "Ada" Wang
# Date: August/21/2020
########################################

library(AugRCT)

##------------- Step 1: load Stan models---------------------------
#specify path of AugRCT package on R server
root.path <- "/SFS/user/ctc/wangxi8/AugRCT_new/"
#options(mc.cores=parallel::detectCores()) #In HPC server, DO NOT run this line of code
rstan_options(auto_write = TRUE) #tell stan not to compile code that has already been compiled
path <- paste0(root.path,"stan/")
# Matching + Fixed power prior. Bayes model 1. logit(theta) =beta0 + lambda*Z
model.1 <- stan_model(paste0(path, "model1.stan"))
# Matching + Fixed power prior. Bayes model 2. theta ~ Beta(k*mu, k*(1-mu))
model.2 <- stan_model(paste0(path, "model2_new.stan"))
# Stratification + random Power Prior
model.3 <- stan_model(paste0(path, "model3.stan"))
# Stratification + fixed Power Prior
model.3fix <- stan_model(paste0(path, "model3fix.stan"))
# Matching + Commensurate prior
model.4 <- stan_model(paste0(path, "model4_new.stan"))
# PS stratification + Commensurate prior
model.5 <- stan_model(paste0(path, "model5_new.stan"))
# IPTW + Commensurate prior
model.6 <- stan_model(paste0(path, "model6_new.stan"))
