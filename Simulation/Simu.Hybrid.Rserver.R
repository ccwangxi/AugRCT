########################################
# Title: One simulation under 32 settings (Hybrid methods)
# Author: Xi "Ada" Wang
# Date: August/21/2020
########################################

library(AugRCT)

##------------- Step 1: load Stan models---------------------------
#specify path of AugRCT package on R server
root.path <- "/SFS/user/ctc/wangxi8/AugRCT_new/"
#options(mc.cores=parallel::detectCores()) #In HPC server, DO NOT run this line of code
rstan_options(auto_write = TRUE) #tell stan not to compile code that has already been compiled
path <- paste0(root.path,"stan\\")
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

##------------- Step 2: Specify Simulation Settings---------------------------
unm.conf.inds <- c(0, # no unmeasured confounding. ps.cov=c("x1","x2","x3","x4")
                   1) # with unmeasured confounding. ps.cov=c("x1","x2","x4")
lambdas <- c(1.185, # with treatment effect
             0) #no treatment effect

#vectors of (x1.mean.CH, x2.p.CH,xi)
xi.1 <- -0.65; xi.2 <- 0.56
x1.mean.CHs=c(40, 40-1.5, # no data source diff
              40, 40-xi.2*10 + 1.5, 40-xi.2*10/2 + 1.5, # with data source diff. thetach > thetaCD
              40, 40-xi.1*10 - 1.5, 40-xi.1*10/2 - 1.5) # with data source diff. thetach < thetaCD
x2.p.CHs=c(0.4, 0.4+0.3, # no data source diff
           0.4, 0.4-0.3, 0.4-0.3,# with data source diff. thetach > thetaCD
           0.4, 0.4+0.3, 0.4+0.3) # with data source diff. thetach < thetaCD
xis <- c(0, 0, # no data source diff
         xi.2, 0, xi.2/2, # with data source diff. thetach > thetaCD
         xi.1, 0, xi.1/2) # with data source diff. thetach < thetaCD
data.source.mx <- data.frame(x1.mean.CHs, x2.p.CHs, xis=xis)

simu_sets <- expand.grid(data.source = 1:nrow(data.source.mx),
                         unm.confs = 1:length(unm.conf.inds),
                         lambdas = lambdas)


#task_id <- 1  #### NEED to specify task_id as a interger when run it on local
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

##------------- Step 3: Run the simulation---------------------------
comps <- list()
#comps.fit <- list()

for(i in 1:32){
  data.source.vec <- data.source.mx[simu_sets$data.source[i],]

  unm.conf <- unm.conf.inds[simu_sets$unm.confs[i]]


  lambda   <- simu_sets$lambdas[i]

  #Fit model
  if (unm.conf ==1) {
    fit <- Simu.Hybrid(seed=task_id, lambda=lambda, set.num = i,
                                xi=data.source.vec$xis, x1.mean.CH=data.source.vec$x1.mean.CHs, x2.p.CH=data.source.vec$x2.p.CHs,
                                ps.covs=c("x2","x3","x4"),
                                n.E = 80, n.CD = 40, n.CH = 300,
                                m2.iter = 40000, m2.adapt_delta=0.99, m4.iter = 20000,m5.iter = 20000)
  }else {
    fit <- Simu.Hybrid(seed=task_id, lambda=lambda, set.num = i,
                                xi=data.source.vec$xis, x1.mean.CH=data.source.vec$x1.mean.CHs, x2.p.CH=data.source.vec$x2.p.CHs,
                                ps.covs=c("x1","x2","x3","x4"),
                                n.E = 80, n.CD = 40, n.CH = 300,
                                m2.iter = 40000, m2.adapt_delta=0.99, m4.iter = 20000,m5.iter = 20000)
  }

  # Save simulated data and simulation results
  comps[[i]] <- fit$comp
  # comps.fit[[i]] <- fit$comp.fits
}

comps.all <- rbindlist(comps)

save(comps.all, file = paste0(task_id, ".Rdata"))
