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
save.freq.res <- function(path.temp,
                          filename,
                          nsim=200,
                          res.path="/SFS/user/ctc/wangxi8/AugRCT_new/Simulation Results/R Data/"){
  # path.temp: path of simulation results in scrach folder
  # filename: name of RData to save the summary of simulation results
  # nsim: number of simulation
  # res.path: path of summary RData we want to save in the main folder on R server
  setwd(path.temp)
  files <- list.files(path = path.temp, pattern = "*.Rdata")
  comps.freq <- list()
  length(files)
  for(iter in 1:nsim){
    load(files[iter]) #comps.freq.one
    comps.freq[[iter]] <- comps.freq.one
  }
  comps.freq <- rbindlist(comps.freq)

  save(comps.freq,file=paste0(res.path,filename))
}

save.freq.res(path.temp="/SFS/scratch/user/freq_RD", nsim=5000, filename="freq_5000.RData")

#--------------------------------------------
# save results from Simu.Hybrid.Rserver.R
save.stan.res <- function(path.temp,
                          filename,
                          nsim=200,
                          res.path="/SFS/user/ctc/wangxi8/AugRCT_new/Simulation Results/R Data/"){
  # path.temp: path of simulation results in scrach folder
  # filename: name of RData to save the summary of simulation results
  # nsim: number of simulation
  # res.path: path of summary RData we want to save in the main folder on R server
  setwd(path.temp)
  files <- list.files(path = path.temp, pattern = "*.Rdata")
  comps <- list()
  length(files)
  for(iter in 1:nsim){
    load(files[iter]) #comps.all
    comps[[iter]] <- comps.all
  }
  comps <- rbindlist(comps)

  save(comps,file=paste0(res.path,filename))
}

save.stan.res(path.temp="/SFS/scratch/user/stan_RD", nsim=5000, filename="stan_5000.RData")
