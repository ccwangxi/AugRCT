#' Implement Hybrid approaches in simulation
#'
#' This function allows you to implement hybrid approaches (PS plus Bayesian) in augmenting current control with historical data in simulation.
#' @param seed a single value, interpreted as an integer, or NULL. seed of random number generator.
#' @param set.num number of simulation setting, range from 1 to 32.
#' @param n.E a positive integer for number of patients in current experimental arm (E). Default to 80.
#' @param n.CD a positive integer for number of patients in current control arm (CD). Default to 40.
#' @param n.CH a positive integer for number of patients in historical control arm (CH). Default to 300.
#' @param lambda log odds of treatment effect. HR=exp(lambda).
#' @param xi log odds of current study effect. = 0 if no drift.
#' @param x1.mean.ECD mean of continuous variable X1 in the current trial. Default to 40
#' @param x1.mean.CH mean of continuous variable X1 in the historical trial. Default to 40
#' @param x2.p.ECD mean of binary variable X2 in the current trial. Default to 0.4
#' @param x2.p.CH mean of binary variable X2 in the historical trial. Default to 0.4
#' @param beta0 intercept in the logistic outcome model. Default to 3.6
#' @param beta1 coefficient of X1 in the logistic outcome model. Default to -0.1
#' @param beta2 coefficient of X2 in the logistic outcome model. Default to -0.5
#' @param beta3 coefficient of X3 in the logistic outcome model. Default to 1
#' @param beta4 coefficient of X4 in the logistic outcome model. Default to -1
#' @param nstrata a positive integer of number of strata in PS stratification. Default to 3.
#' @param trt a character of variables name for the binary outcome in PS model. Defaul to "D"
#' @param ps.covs a list of characters for variable na,es included in PS model. Default to \code{c("x1","x2","x3","x4")}
#' @param m2.iter a positive integer of iteration number in STAN model 2 (Matching + Fixed power prior. Beta model). Default to 40000
#' @param m2.adapt_delta a real number between 0 and 1 to control the step size in STAN model 2. Default to 0.99
#' @param m2.max_treedepth a positive integer of tree depth of STAN model 2. Default to 12
#' @param m4.iter a positive integer of iteration number in STAN model 4 (Matching + Commensurate prior). Default to 20000
#' @param m4.adapt_delta a real number between 0 and 1 to control the step size in STAN model 4. Default to 0.99
#' @param m5.iter a positive integer of iteration number in STAN model 5 (PS stratification + Commensurate prior). Default to 20000
#' @param m5.adapt_delta a real number between 0 and 1 to control the step size in STAN model 5. Default to 0.99
#' @param seed a single value, interpreted as an integer, or NULL. seed of random number. Default to 2020
#' @return mydata, a data frame of simulated dataset
#' @return comp, a data frame of the point estimate and percentiles of delta (response rate difference between experiment and augmented control group) of all the hybrid methods.
#' @keywords AugRCT
#' @import dplyr tidyr MatchIt rstan
#' @details There might be warnings on the divergence of the STAN model fit. But the chains of posterior samples of the parameter of interest mixes well by checking the trace plot.
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu
#' @export
#'
Simu.Hybrid <- function(seed, set.num=0,
                        n.E = 80, n.CD = 40, n.CH = 300,
                        lambda=2,xi=0,
                        x1.mean.ECD=40,x1.mean.CH=40,
                        x2.p.ECD=0.4,x2.p.CH=0.4,
                        beta0=3.6,beta1=-0.1,beta2=-0.5,beta3=1,beta4=-1,
                        nstrata = 3,
                        m2.iter = 40000, m2.adapt_delta=0.99, m2.max_treedepth = 12,
                        m4.iter = 20000, m4.adapt_delta=0.99,
                        m5.iter = 20000, m5.adapt_delta=0.99,
                        trt="D",ps.covs=c("x1","x2","x3","x4")){
  mydata <- generate.data(seed=seed,n.E=n.E,n.CD=n.CD,n.CH=n.CH,lambda=lambda,xi=xi,
                          x1.mean.ECD=x1.mean.ECD,x1.mean.CH=x1.mean.CH,
                          x2.p.ECD=x2.p.ECD,x2.p.CH=x2.p.CH,
                          beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4)
  # get sample.theta.E, sample.theta.CD, sample.theta.CH.
  mydata.ps <- add.ps(data=mydata,ps.covs=ps.covs)
  mydata <- mydata.ps$data.ps

  #pair matching
  PS.distance <- match_on(D ~ est.PS, method="euclidean",data=mydata)
  pairmatch <- pairmatch(PS.distance, data=mydata)
  #Generate Matched dataset. randomly select n.E - n.CD from n.E matched subjects
  mydata1 <- data.frame(mydata, matches = pairmatch, check.rows=T)
  pm.mydata <- mydata1 %>% dplyr::filter(group=="CD" | (!is.na(matches))) %>%
    dplyr::mutate(keep=ifelse(group %in% c("E","CD"),1,
                              sample(c(rep(0,n.CD*2),rep(1,n.E - n.CD)),replace=FALSE))) %>%
    dplyr::filter(keep==1)

  ds <- pm.mydata %>% arrange(CH.ind) #order: E CD CH
  data <- list(NECD=n.CD+n.E,N=n.E*2,
               a0=ds$est.PS, Y=ds$Y, Z=ds$Z)
  fit.1 <- rweSTAN(lst.data = data, model = model.1)

  data <- list(NECD = n.CD+n.E,N = n.E*2,NC = n.E,NCH = n.E-n.CD,
               a0 = ds$est.PS, Y = ds$Y, Z = ds$Z)
  fit.2 <- rweSTAN(lst.data = data, model = model.2,iter = m2.iter,control = list(adapt_delta=m2.adapt_delta,max_treedepth = m2.max_treedepth))

  # IPTW + fixed power prior with (stabilized IPTW (ATT) weight)
  prop <- mydata %>% dplyr::group_by(D) %>% dplyr::summarise(mean = mean(Y))
  suppressMessages({mydata.iptw <- mydata %>%
    dplyr::mutate(weight.ns = ifelse(mydata$group %in% c("E","CD"),1, est.PS/(1-est.PS)))})# weight.ns: non-stabilized weight
  mydata.iptw$weight[mydata.iptw$D==0] <- mydata.iptw$weight.ns[mydata.iptw$D==0]/mean(mydata.iptw$weight.ns[mydata.iptw$D==0])
  if (max(mydata.iptw$weight) > 10) {print("Stabilized IPTW weight > 10")}
  ds <- mydata.iptw %>% arrange(CH.ind) #order: E CD CH
  data <- list(NECD=n.CD+n.E,N=dim(ds)[1],
               a0=ds$weight*ds$est.PS, Y=ds$Y, Z=ds$Z)
  fit.6 <- rweSTAN(lst.data = data, model = model.1)

  # IPTW + commensurate prior
  ds <- ds %>% dplyr::arrange(CH.ind,Z) # order of group: CD, E CH
  data <- list(NE = n.E,NCD = n.CD,N = dim(ds)[1],Y = ds$Y,prior_alpha = 1,prior_beta = 1,a0=ds$weight)
  fit.7 <- rweSTAN(lst.data = data, model = model.6,iter = m5.iter,control = list(adapt_delta=m5.adapt_delta))

  fit.trim.1 <- function(data){
    #this function is used to generate the data for IPTW trimming analysis with weight
    #1. find delta.trim
    data <- data %>%
      dplyr::mutate(weight = ifelse(data$group %in% c("E","CD"),1, est.PS/(1-est.PS)))
    find.trim.fn <- function(x){
      #this function is used to find the delta.trim for IPTW trimming of data
      abs(mean(data$Y[data$group == "CD"]) -
            sum(((data$Y)*(data$weight))[data$group == "CH" & (data$weight>x)])/sum((data$weight)[data$group == "CH" & (data$weight>x)]))
    }
    # get the value of lambda
    min <- min(data$weight[data$group == "CH"])+1e-3; max <- max(data$weight[data$group == "CH"])
    delta.trim <- optimize(find.trim.fn,maximum = FALSE,tol=1e-7,interval=c(min,max))$minimum
    #2. get trimmed data
    data.trim <- data %>% dplyr::mutate(keep=ifelse(group %in% c("CD","E"),1,
                                                    ifelse(weight > delta.trim,1,0))) %>% dplyr::filter(keep==1)
    data.trim$weight.ns <- data.trim$weight
    data.trim$weight[data.trim$D==0] <- data.trim$weight.ns[data.trim$D==0]/mean(data.trim$weight.ns[data.trim$D==0])
    return(data.trim)
  }
  # IPTW trimming + fixed power prior
  ds <- fit.trim.1(mydata) %>% arrange(CH.ind) #order: E CD CH
  data <- list(NECD=n.CD+n.E,N=dim(ds)[1],
               a0=ds$weight*ds$est.PS, Y=ds$Y, Z=ds$Z)
  fit.8 <- rweSTAN(lst.data = data, model = model.1)

  # Matching + Commensurate prior
  ds <- pm.mydata %>% dplyr::arrange(CH.ind,Z) # order of group: CD, E CH
  data <- list(NE = n.E,NCD = n.CD,NCH = n.E-n.CD,Y = ds$Y,prior_alpha = 1,prior_beta = 1)
  fit.4 <- rweSTAN(lst.data = data, model = model.4,iter = m4.iter,control = list(adapt_delta=m4.adapt_delta))

  #IPTW + trimming + commensurate prior
  ds <- fit.trim.1(mydata) %>% dplyr::arrange(CH.ind,Z) # order of group: CD, E CH
  data <- list(NE = n.E,NCD = n.CD,N = dim(ds)[1],Y = ds$Y,prior_alpha = 1,prior_beta = 1,a0=ds$weight)
  fit.9 <- rweSTAN(lst.data = data, model = model.6,iter = m5.iter,control = list(adapt_delta=m5.adapt_delta))

  data.withps <- rwePS(mydata, ps.fml = NULL,v.grp = trt, #treatment variable
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
  # PS Stratification + commensurate prior
  rwe.data[["prior_alpha"]]=1;rwe.data[["prior_beta"]]=1
  fit.5 <- rweSTAN(lst.data = rwe.data,
                   model = model.5,iter = m5.iter,control = list(adapt_delta=m5.adapt_delta));

  fits <- list(fit.1, fit.2, fit.4, fit.31, fit.32, fit.5, fit.6, fit.7, fit.8, fit.9)

  labels <- c("matching + power prior (Bayes model 1)","matching + power prior (Bayes model 2)",
              "Matching + Commensurate prior","Stratification + random Power Prior",
              "Stratification + fixed Power Prior","Stratification + commensurate prior",
              "IPTW + fixed power prior (logit model)",
              "IPTW + commensurate prior",
              "IPTW trim + fixed power prior (logit model)",
              "IPTW trim + commensurate prior")#"nearest matching + power prior (Bayes model 1)",
  comp <- comp.fits(fits, labels)
  comp <- comp %>% dplyr::mutate(sim = seed, set.num = set.num)
  return(list(mydata=mydata,comp=comp))
}
