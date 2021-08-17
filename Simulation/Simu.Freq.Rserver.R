########################################
# Title: One simulation under 32 settings (Frequentist methods)
# Author: Xi "Ada" Wang
# Date: August/21/2020
########################################

library(AugRCT)

##------------- Step 1: Specify Simulation Settings---------------------------
unm.conf.inds <- c(0, # no unmeasured confounding. ps.cov=c("x1","x2","x3","x4")
                   1) # with unmeasured confounding. ps.cov=c("x1","x2","x4")
lambdas <- c(1.185, # with treatment effect
             0) #no treatment effect

#vectors of (x1.mean.CH, x2.p.CH,xi)
xi.1 <- -0.65; xi.2 <- 0.56   #### CHANGE this line for different drift size. The value here corresponds to a drift size of 0.1
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

##------------- Step 2: Run the simulation---------------------------
comps.freq <- list()

for(i in 1:32){

  data.source.vec <- data.source.mx[simu_sets$data.source[i],]

  unm.conf <- unm.conf.inds[simu_sets$unm.confs[i]]


  lambda   <- simu_sets$lambdas[i]

  #Fit model
  if (unm.conf ==1) {
    fit <- Simu.Freq(seed=task_id, lambda=lambda, set.num = i,
                                xi=data.source.vec$xis, x1.mean.CH=data.source.vec$x1.mean.CHs, x2.p.CH=data.source.vec$x2.p.CHs,
                                ps.covs=c("x2","x3","x4"),
                                n.E = 80, n.CD = 40, n.CH = 300,
                                m2.iter = 40000, m2.adapt_delta=0.99, m4.iter = 20000,m5.iter = 20000)
  }else {
    fit <- Simu.Freq(seed=task_id, lambda=lambda, set.num = i,
                                xi=data.source.vec$xis, x1.mean.CH=data.source.vec$x1.mean.CHs, x2.p.CH=data.source.vec$x2.p.CHs,
                                ps.covs=c("x1","x2","x3","x4"),
                                n.E = 80, n.CD = 40, n.CH = 300,
                                m2.iter = 40000, m2.adapt_delta=0.99, m4.iter = 20000,m5.iter = 20000)
  }

  # Save simulated data and simulation results
  comps.freq[[i]] <- fit$comp.freq

}

data <- fit$mydata
comps.freq.one <- rbindlist(comps.freq)

save(data, comps.freq.one, file = paste0(task_id, ".Rdata"))
