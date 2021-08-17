########################################
# Title: Compariosn of Subgroup Methods in Simulation for Augmented control
# Author: Xi "Ada" Wang
# Date: Aug/21/2020
########################################

library(knitr)
library(markdown)
library(rmarkdown)
library(stringr)

root.path <- "C:\\Users\\wangxi8\\Documents\\AugRCT\\"
library(AugRCT)

#-------------------step 1: Prepare Data-------------------------------------
## Simplify Simulation RD (All, No DS, DS larger, DS smaller)
load(file = paste0(root.path,"Simulation results\\R Data\\simu_true_theta1.RData")) #unm.conf.inds, data.source.mx, simu.true
load(file = paste0(root.path,"Simulation results\\R Data\\freq_RD_500_0820.RData"))
unique(comps.freq$label)
comps.freq <- comps.freq %>%
  dplyr::filter(! label %in% c("standard"))

load(file = paste0(root.path,"Simulation results\\R Data\\stan_RD_500_0803.RData"))
unique(comps$fit)
comps <- comps %>% dplyr::filter(! fit %in% c("matching + power prior (Bayes model 1)","Stratification + fixed Power Prior",
                                              "IPTW trim + fixed power prior (logit model)", "IPTW trim + commensurate prior"))

comps.freq$set.num <- as.character(comps.freq$set.num)
comps$set.num <- as.character(comps$set.num)
comps <- merge(comps, simu.true, by="set.num")
comps.freq <- merge(comps.freq, simu.true, by="set.num")


# generate table
comps <- comps %>% dplyr::mutate(point.est = mean, lower.CI = X2.5., upper.CI = X97.5.,label = fit) %>%
  dplyr::select(point.est,lower.CI,upper.CI,label,sim,
                set.num, data.source, unm.confs, lambdas, true.thetaE, true.thetaCD,true.thetaCH, data.source.diff,true.delta)

method.level <- c("No borrowing", "Pool",
                  #"Matching for order","PS stratification for order", "IPTW for order", #add these three levels for order within group
                  "Matching", "matching + power prior (Bayes model 2)","Matching + Commensurate prior",
                  "PS stratification (trim)","Stratification + random Power Prior", "Stratification + commensurate prior",
                  "IPTW", "IPTW trim (org)","IPTW + fixed power prior (logit model)","IPTW + commensurate prior")
method.label <- c("No borrowing", "Pooling",
                  # "PS Matching","PS stratification", "IPTW",
                  "PS Matching", "Power Prior","Commensurate Prior",
                  "PS Stratification","Power prior", "Commensurate prior",
                  "IPTW","Trimming","power prior", "commensurate prior")
subgroup.label <- c(rep("Frequentist",2),rep("PS Matching Only/+ Bayesian \n Information Borrowing", 3),
                    rep("PS Stratification Only/+ Bayesian \n Information Borrowing", 3),
                    rep("IPTW Only/+Trimming/+Bayesian \n Information Borrowing", 4))
dt.all <- rbind(comps, comps.freq)


cols = c("goldenrod1","gold3","olivedrab3","seagreen","springgreen3","#56B4E9","#0072B2","deepskyblue3")

dts <- list(dt.all, dt.all %>% dplyr::filter(data.source %in% c(1:2)),
            dt.all %>% dplyr::filter(data.source %in% c(3:5)),dt.all %>% dplyr::filter(data.source %in% c(6:8)))
cols.list <- list(cols, cols[1:2], cols[3:5], cols[6:8])
method.levels <- list(method.level)
method.labels <- list(method.label)
subgroup.labels <- list(subgroup.label)
Dir.names <- list("RD Simplify","RD Simplify No DS","RD Simplify Larger","RD Simplify Smaller")
unique(dt.all$label)
n.sim = 500

#-------------------step 2: knitr loop-------------------------------------
for (i in 1:4){
  Dir.name <- Dir.names[[i]]
  dt.all <- dts[[i]]
  method.level <- method.levels[[1]]
  method.label <- method.labels[[1]]
  subgroup.label <- subgroup.labels[[1]]
  cols <- cols.list[[i]]
  rmarkdown::render(paste0(root.path, "Simulation Results\\RMarkdown All Comparison.Rmd"),
                    output_file = paste(Dir.name, "_", format(Sys.Date(),format="%m%d%Y"), ".html",sep=""),
                    output_dir = paste0(root.path, "Simulation Results\\"))
}
