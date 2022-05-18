########################################
# Title: Summary simulation results (Frequentist/Hybrid methods)
# Author: Xi "Ada" Wang
# Date: August/21/2020
########################################
# Perform Sumarry Analysis in the Scrach Folder and save the publication ready simulation tables
# simulation_result/folder.

library(data.table)
library(dplyr)

#--------------------------------------------
# save results from Simu.Freq.Rserver.R
save.freq.res(path.temp="/SFS/scratch/user/freq_RD", nsim=5000, filename="freq_5000.RData")

#--------------------------------------------
# save results from Simu.Hybrid.Rserver.R
save.stan.res(path.temp="/SFS/scratch/user/stan_RD", nsim=5000, filename="stan_5000.RData")
