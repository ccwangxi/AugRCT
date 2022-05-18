########################################
# Title: Summarize Simulation Results
# Author: Xi "Ada" Wang
# Date: August/21/2020
########################################
library(knitr)
library(markdown)
library(rmarkdown)
library(stringr)
# load R functions
root.path <- "H:\\research\\AugRCT\\"
code.path <- paste0(root.path,"R\\")
source(paste0(code.path,"RFunctions.R"))
cols = c("goldenrod1","gold3","olivedrab3","seagreen","springgreen3","#56B4E9","#0072B2","deepskyblue3")

# Step 1: Prepare Clean data ---------
load(file = paste0(root.path,"Simulation results\\R Data\\simu_true_theta1.RData")) #see Save.Simulation.Settings.R file
load(file = paste0(root.path,"Simulation results\\R Data\\freq_5000.RData"))
n.sim = 5000
comps.freq$set.num <- as.character(comps.freq$set.num)
comps.freq <- merge(comps.freq, simu.true, by="set.num")
load(file = paste0(root.path,"Simulation results\\R Data\\stan_5000.RData"))
comps$set.num <- as.character(comps$set.num)
comps <- merge(comps, simu.true, by="set.num")
addcov.setting.ind = 0

# Step 2: Calculate metrics for simulations ---------
RD.bound <- 0
if (addcov.setting.ind == 1){
  Drift.labels <- c("1. 0 (No CDD/Time Trend)","2. 0 (CDD)","3. 0.2 (Time trend)", "4. 0.2 (CDD)", "5. 0.2 (Time trend + CDD)", "6. -0.2 (Time trend)", "7. -0.2 (CDD)", "8. -0.2 (Time trend + CDD)")
}else{
  Drift.labels <- c("1. 0 (No CDD/Time Trend)","2. 0 (CDD)","3. 0.1 (Time trend)", "4. 0.1 (CDD)", "5. 0.1 (Time trend + CDD)", "6. -0.1 (Time trend)", "7. -0.1 (CDD)", "8. -0.1 (Time trend + CDD)")
}

comps.all <- dt.all %>% dplyr::mutate(true.par = true.delta) %>%
  dplyr::mutate(cov.ind = 1*((true.par >= lower.CI) & (true.par <= upper.CI)),
                Power.ind = 1*(lower.CI > RD.bound),
                Bias.ind = point.est - true.par,
                CI.width = upper.CI - lower.CI)

# summ include simulation results
summ <- comps.all %>% dplyr::group_by(set.num, label, data.source, unm.confs, lambdas) %>%
  dplyr::summarise(nsim = n(),
                   Bias = mean(Bias.ind),
                   Coverage = 100*mean(cov.ind),
                   Emp.SD = sd(point.est),
                   Emp.MSE = sqrt(sum(Bias.ind^2)/(nsim-1)),
                   Emp.MSE2 = sum(Bias.ind^2)/(nsim-1),
                   CI.width.mean = mean(CI.width),
                   Power = mean(Power.ind))

