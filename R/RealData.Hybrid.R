#' Implement Hybrid approaches in real data application
#'
#' This function allows you to implement hybrid approaches (PS plus Bayesian) in augmenting current control with historical data in real data application.
#' @param realdata the name of a data frame containing all needed variables, including D (indicator of current trial), Z (indicator of treatment),
#' Y (binary outcome), a list of baseline covariates included in PS model
#' @param ps.covs a list of characters for variable na,es included in PS model
#' @param trt a character of variables name for the binary outcome in PS model. Defaul to "D"
#' @param nstrata a positive integer of number of strata in PS stratification. Default to 5.
#' @param m2.iter a positive integer of iteration number in STAN model 2 (Matching + Fixed power prior. Beta model). Default to 60000
#' @param m2.adapt_delta a real number between 0 and 1 to control the step size in STAN model 2. Default to 0.999
#' @param m2.max_treedepth a positive integer of tree depth of STAN model 2. Default to 12
#' @param m4.iter a positive integer of iteration number in STAN model 4 (Matching + Commensurate prior). Default to 20000
#' @param m4.adapt_delta a real number between 0 and 1 to control the step size in STAN model 4. Default to 0.99
#' @param m5.iter a positive integer of iteration number in STAN model 5 (PS stratification + Commensurate prior). Default to 20000
#' @param m5.adapt_delta a real number between 0 and 1 to control the step size in STAN model 5. Default to 0.99
#' @param seed a single value, interpreted as an integer, or NULL. seed of random number. Default to 2020
#' @return comp, a data frame of the point estimate and percentiles of delta (response rate difference between experiment and augmented control group) of all the hybrid methods.
#' @return comp.fits, the list of the model fits for all the hybrid methods.
#' @keywords AugRCT
#' @import dplyr tidyr MatchIt rstan
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu
#' @export


# load functions under "stan" folder
RealData.Hybrid <- function(realdata,
                            ps.covs, trt="D", nstrata = 5,
                            m2.iter = 60000, m2.adapt_delta=0.999, m2.max_treedepth = 12,
                            m4.iter = 20000, m4.adapt_delta=0.99,
                            m5.iter = 20000, m5.adapt_delta=0.99,
                            seed=2020){
  set.seed(seed)
  n.E = sum(realdata$group=="E"); n.CD = sum(realdata$group=="CD"); n.CH = sum(realdata$group=="CH")

  # add estimated PS
  realdata.ps <- Add.PS(data=realdata, trt=trt, ps.covs=ps.covs)
  realdata$est.PS <- realdata.ps$data.ps$est.PS

  #pair matching
  PS.distance <- match_on(D ~ est.PS, method="euclidean",data=realdata)
  pairmatch <- pairmatch(PS.distance, data=realdata)
  #Generate Matched dataset. randomly select n.E - n.CD from n.E matched subjects
  realdata1 <- data.frame(realdata, matches = pairmatch, check.rows=T)
  if (n.E + n.CD < n.CH) {
    keep.num <- n.E - n.CD; drop.num <- n.E + n.CD - (n.E - n.CD)
  }else {
    keep.num <- n.E - n.CD; drop.num <- n.CH - (n.E - n.CD)
  }
  pm.realdata.check <- realdata1 %>%
    dplyr::mutate(keep=ifelse(group %in% c("E","CD"),1,
                              ifelse(is.na(matches), 0, NA)))
  random.list <- sample(c(rep(0,drop.num),rep(1,keep.num)),replace=FALSE)
  pm.realdata.check$keep[is.na(pm.realdata.check$keep)] <- random.list
  pm.realdata <- pm.realdata.check %>%
    dplyr::filter(keep==1)

  ds <- pm.realdata %>% arrange(desc(D)) #order: E CD CH
  data <- list(NECD=n.CD+n.E,N=n.E*2,
               a0=ds$est.PS, Y=ds$Y, Z=ds$Z)
  fit.1 <- rweSTAN(lst.data = data, model = model.1)

  data <- list(NECD = n.CD+n.E,N = n.E*2,NC = n.E,NCH = n.E-n.CD,
               a0 = ds$est.PS, Y = ds$Y, Z = ds$Z)
  fit.2 <- rweSTAN(lst.data = data, model = model.2,iter = m2.iter,control = list(adapt_delta=m2.adapt_delta,max_treedepth = m2.max_treedepth))

  ds <- pm.realdata %>% dplyr::arrange(desc(D),Z) # order of group: CD, E CH
  data <- list(NE = n.E,NCD = n.CD,NCH = n.E-n.CD,Y = ds$Y,prior_alpha = 1,prior_beta = 1)
  fit.4 <- rweSTAN(lst.data = data, model = model.4,iter = m4.iter,control = list(adapt_delta=m4.adapt_delta))

  # IPTW + fixed power prior with (stabilized IPTW (ATT) weight)
  prop <- realdata %>% dplyr::group_by(D) %>% dplyr::summarise(mean = mean(Y)); prop.D = prop$mean[prop$D==1]
  suppressMessages({realdata.iptw <- realdata %>%
    dplyr::mutate(weight = ifelse(realdata$group %in% c("E","CD"),1, (1-prop.D)*est.PS/(prop.D*(1-est.PS))))})
  ds <- realdata.iptw %>% arrange(desc(D)) #order: E CD CH
  data <- list(NECD=n.CD+n.E,N=dim(ds)[1],
               a0=ds$weight*ds$est.PS, Y=ds$Y, Z=ds$Z)
  fit.6 <- rweSTAN(lst.data = data, model = model.1)

  # IPTW + commensurate prior
  ds <- ds %>% dplyr::arrange(desc(D),Z) # order of group: CD, E CH
  data <- list(NE = n.E,NCD = n.CD,N = dim(ds)[1],Y = ds$Y,prior_alpha = 1,prior_beta = 1,a0=ds$weight)
  fit.7 <- rweSTAN(lst.data = data, model = model.6,
                   iter = m5.iter,control = list(adapt_delta=m5.adapt_delta))

  data.withps <- rwePS(realdata, ps.fml = NULL,v.grp = trt, #treatment variable
                       v.covs = ps.covs, #list of covariates
                       d1.grp = 1, #level of D
                       nstrata = nstrata)
  psdist <- rwePSDist(data.withps,
                      n.bins = 10, # number of bins to categorize distribution of PS
                      type = "ovl")
  rwe.data <- rwePsPreprocess(data = data.withps$data,
                              A = n.E-n.CD, #target number of subjects to be borrowed
                              RS = psdist[psdist$Strata %in% 1:nstrata,"Dist"],
                              v.outcome = "Y")
  # PS stratification + random Power Prior
  fit.31 <- rweSTAN(lst.data = rwe.data, model = model.3);
  # PS stratification + fixed Power Prior
  fit.32 <- rweSTAN(lst.data = rwe.data, model = model.3fix);
  rwe.data[["prior_alpha"]]=1;rwe.data[["prior_beta"]]=1
  # PS Stratification + commensurate prior
  fit.5 <- rweSTAN(lst.data = rwe.data, model = model.5,iter = m5.iter,control = list(adapt_delta=m5.adapt_delta));

  fits <- list(fit.1, fit.2, fit.4,
               fit.31, fit.32, fit.5,
               fit.6, fit.7)

  labels <- c("matching + power prior (Bayes model 1)","matching + power prior (Bayes model 2)","Matching + Commensurate prior",
              "Stratification + random Power Prior","Stratification + fixed Power Prior","Stratification + commensurate prior",
              "IPTW + fixed power prior (logit model)", "IPTW + commensurate prior")
  comp <- Extract.Hybrid(fits, labels)
  return(list(comp=comp, comp.fits=fits))
}
