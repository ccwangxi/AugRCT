########################################
# Title: Save Simulation Settings
# Author: Xi "Ada" Wang
# Date: August/21/2020
########################################

##------------- Step 1: Calculate Parameters for Simulation Settings---------------------------
lambda <- find.lambda(n.sim=1000,power=0.8)
lambda #1.185
xi.1 <- find.xi(thetaCH=0.317-0.1) #-0.6459937
xi.1
xi.2 <- find.xi(thetaCH=0.317+0.1) #0.5588587
xi.2

xi.1 <- find.xi(thetaCH=0.317-0.2) #-1.529172
xi.1
xi.2 <- find.xi(thetaCH=0.317+0.2) #1.084294
xi.2

get.simu.setting <- function(lambda,x1.mean.CH=40, x2.p.CH=0.4, xi=0,
                             n.E=80, n.CD=40, n.CH=300,
                             x1.mean.ECD=40,
                             x2.p.ECD=0.4,
                             beta0=3.6, beta1=-0.1, beta2=-0.5, beta3=1, beta4=-1){
  # this function is used to generate true.delta for simulated datasets
  mydata <- generate.data(seed=2020, n.E=n.E*10000, n.CD=n.CD*10000, n.CH=n.CH*10000,
                          lambda=lambda, xi=xi,
                          x1.mean.ECD=x1.mean.ECD, x1.mean.CH=x1.mean.CH,
                          x2.p.ECD=x2.p.ECD, x2.p.CH=x2.p.CH,
                          beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4)
  trues <- mydata %>% dplyr::group_by(group) %>% dplyr::summarise(mean=mean(Y))
  true.thetaE = trues$mean[trues$group=="E"]
  true.thetaCD = trues$mean[trues$group=="CD"]
  true.thetaCH = trues$mean[trues$group=="CH"]
  true.delta = true.thetaE - true.thetaCD
  trues.out <- round(data.frame(true.thetaE, true.thetaCD, true.thetaCH, true.delta),3)
  trues.out
}

##------------- Step 2: Specify Simulation Settings---------------------------
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

##------------- Step 2.5:Calculate true.delta for Simulation Settings---------------------------
trues.delta <- list()
for(i in 1:nrow(simu_sets)){
  data.source.vec <- data.source.mx[simu_sets$data.source[i],]

  unm.conf <- unm.conf.inds[simu_sets$unm.confs[i]]
  if (unm.conf ==1) {
    ps.cov <- ps.cov <- c("x2","x3","x4")
  }else {ps.cov <- c("x1","x2","x3","x4")}

  lambda   <- simu_sets$lambdas[i]

  trues <- get.simu.setting(x1.mean.CH=data.source.vec$x1.mean.CHs,
                            x2.p.CH=data.source.vec$x2.p.CHs,
                            xi=data.source.vec$xis,
                            lambda=lambda,
                            addcov.setting.ind = 1)
  trues.delta[[i]] <- trues
  set.num = i
}

simu.true0 <- data.frame(
  simu_sets,
  do.call(rbind,trues.delta))

simu.true <- simu.true0 %>%
  dplyr::mutate(data.source.diff = true.thetaCD - true.thetaCH,
                set.num = rownames(simu.true)) %>%
  dplyr::select(set.num, data.source, unm.confs, lambdas, true.thetaE, true.thetaCD, true.thetaCH, data.source.diff, true.delta)

# save simulation settings information
root.path <- "H:\\research\\AugRCT\\"
save(simu.true, file = paste0(root.path, "Simulation results\\R Data\\simu_true_theta1.RData"))
