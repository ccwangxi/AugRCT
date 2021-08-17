#' Implement Frequentist (naive and PS only) approaches in simulation
#'
#' This function allows you to implement all the frequentist approaches (naive and PS only) in augmenting current control with historical data in simulation.
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
#' @param bootstrap.seed seed for bootstrap CI of IPTW approach. Default to 1234
#' @param bootstrap.rep replication number of Bootstrap samples. Default to 5000
#' @return mydata, a data frame of simulated dataset
#' @return comp.freq, a data frame of the point estimate and percentiles of delta (response rate difference between experiment and augmented control group) of all the frequentist methods.
#' @keywords AugRCT
#' @import dplyr tidyr MatchIt rstan
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu
#' @export

Simu.Freq <- function(seed, set.num=0,
                      n.E = 80, n.CD = 40, n.CH = 300,
                      lambda=2,xi=0,
                      x1.mean.ECD=40,x1.mean.CH=40,
                      x2.p.ECD=0.4,x2.p.CH=0.4,
                      beta0=3.6,beta1=-0.1,beta2=-0.5,beta3=1,beta4=-1,
                      nstrata = 3,
                      trt="D",ps.covs=c("x1","x2","x3","x4"),
                      bootstrap.seed=1234, bootstrap.rep=5000){
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
  #table(!is.na(mydata1$matches),mydata1$group)
  pm.mydata <- mydata1 %>% dplyr::filter(group=="CD" | (!is.na(matches))) %>%
    dplyr::mutate(keep=ifelse(group %in% c("E","CD"),1,
                              sample(c(rep(0,n.CD*2),rep(1,n.E - n.CD)),replace=FALSE))) %>%
    dplyr::filter(keep==1)

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

  #------------------------------------------------
  # Frequentist approaches
  #PS Matching + logistic outcome model
  outcome1 <- get.mn.ci(data = pm.mydata, label = "Matching")

  # IPTW (ATT)  + logistic outcome model (MSM marginal structure model. no truncation of weights for now)
  prop <- mydata %>% dplyr::group_by(D) %>% dplyr::summarise(mean = mean(Y)); prop.D = prop$mean[prop$D==1]
  suppressMessages({mydata.iptw <- mydata %>%
    dplyr::mutate(weight.ns = ifelse(mydata$group %in% c("E","CD"),1, est.PS/(1-est.PS)))})# weight.ns: non-stabilized weight
  mydata.iptw$weight[mydata.iptw$D==0] <- mydata.iptw$weight.ns[mydata.iptw$D==0]/mean(mydata.iptw$weight.ns[mydata.iptw$D==0])
  if (max(mydata.iptw$weight) > 10) {print("Stabilized IPTW weight > 10")}

  data <- mydata.iptw
  ps.fml <- as.formula(paste("D", "~", paste(c(ps.covs), collapse="+"), sep=""))
  ps.fit <- glm(ps.fml, data = data, family = "binomial")
  Prop.Z = mean(data$Z)
  data <- data %>%
    dplyr::mutate(delta11 = Z*D*Y/Prop.Z ,
                  delta01 = (1-Z)*D*Y/(1-Prop.Z),
                  delta00 = weight*(1-Z)*(1-D)*Y/(1-Prop.Z))
  delta = mean(data$delta11 - (data$delta01 + data$delta00))
  # bootstrap variance, percentile
  set.seed(1234)
  B <- bootstrap.rep # bootstrap times
  N.tot <- dim(data)[1]
  Boot.rows <- matrix(sample(1:N.tot, size=B*N.tot, replace=TRUE),
                      nrow=N.tot, ncol=B)
  Boot.deltas <- apply(Boot.rows, MARGIN = 2, FUN=function(x){
    data <- mydata.iptw[x,]
    ps.fml <- as.formula(paste("D", "~", paste(c(ps.covs), collapse="+"), sep=""))
    ps.fit <- glm(ps.fml, data = data, family = "binomial")
    Prop.Z = mean(data$Z)
    data <- data %>%
      dplyr::mutate(delta11 = Z*D*Y/Prop.Z ,
                    delta01 = (1-Z)*D*Y/(1-Prop.Z),
                    delta00 = weight*(1-Z)*(1-D)*Y/(1-Prop.Z))
    delta = mean(data$delta11 - (data$delta01 + data$delta00))
    return(delta)
  })
  res.tmp <- data.frame(point.est=delta, lower.CI=quantile(Boot.deltas, probs=0.025), upper.CI=quantile(Boot.deltas, probs=0.975))
  outcome2 <- cbind(res.tmp, label="IPTW")

  # logistic outcome model (Pool, full-borrowing)
  outcome3 <- get.mn.ci(data = mydata, label = "Pool")
  # logistic outcome model (no borrowing)
  outcome4 <- get.mn.ci(data = mydata %>% dplyr::filter(group != "CH"), label = "No borrowing")
  # logistic outcome model (standard setting. RCT 1:1. 80 CD and 80 E)
  add.CD <- generate.data(seed = seed+1,n.E = n.E,n.CD = n.CD,n.CH = n.CH, #use different seed to generate additional n.CD data
                          lambda = lambda,xi = xi,
                          x1.mean.ECD = x1.mean.ECD,x1.mean.CH = x1.mean.CH,
                          x2.p.ECD = x2.p.ECD,x2.p.CH = x2.p.CH,
                          beta0 = beta0,beta1 = beta1,beta2 = beta2,beta3 = beta3,beta4 = beta4)
  standard.data <- rbind(mydata %>% dplyr::filter(group %in% c("E","CD")) %>% dplyr::select(-est.PS),
                         add.CD %>% dplyr::filter(group == "CD"))
  outcome6 <- get.mn.ci(data = standard.data, label = "standard")

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

  # IPTW trimming. original scale
  data <- fit.trim.1(data=mydata)
  ps.fml <- as.formula(paste("D", "~", paste(c(ps.covs), collapse="+"), sep=""))
  ps.fit <- glm(ps.fml, data = data, family = "binomial")
  Prop.Z = mean(data$Z)
  data <- data %>%
    dplyr::mutate(delta11 = Z*D*Y/Prop.Z ,
                  delta01 = (1-Z)*D*Y/(1-Prop.Z),
                  delta00 = weight*(1-Z)*(1-D)*Y/(1-Prop.Z))
  delta = mean(data$delta11 - (data$delta01 + data$delta00))
  # bootstrap variance, percentile
  set.seed(bootstrap.seed)
  B <- bootstrap.rep # bootstrap times
  N.tot <- dim(data)[1]
  Boot.rows <- matrix(sample(1:N.tot, size=B*N.tot, replace=TRUE),
                      nrow=N.tot, ncol=B)
  Boot.deltas <- apply(Boot.rows, MARGIN = 2, FUN=function(x){
    data <- mydata.iptw[x,]
    ps.fml <- as.formula(paste("D", "~", paste(c(ps.covs), collapse="+"), sep=""))
    ps.fit <- glm(ps.fml, data = data, family = "binomial")
    Prop.Z = mean(data$Z)
    data <- data %>%
      dplyr::mutate(delta11 = Z*D*Y/Prop.Z ,
                    delta01 = (1-Z)*D*Y/(1-Prop.Z),
                    delta00 = weight*(1-Z)*(1-D)*Y/(1-Prop.Z))
    delta = mean(data$delta11 - (data$delta01 + data$delta00))
    return(delta)
  })
  res.tmp <- data.frame(point.est=delta, lower.CI=quantile(Boot.deltas, probs=0.025), upper.CI=quantile(Boot.deltas, probs=0.975))
  outcome12 <- cbind(res.tmp, label="IPTW trim (org)")

  # PS stratification
  data <- mydata
  all.ps  <- data$est.PS;
  D1.ps   <- all.ps[which(1 == data[[trt]])]; # CD + E
  cuts    <- quantile(D1.ps, seq(0, 1,length=nstrata+1));
  cuts[1] <- cuts[1] - 0.0001;
  strata.trim <- rep(NA, length(all.ps));
  for (i in 2:length(cuts)) {
    inx      <- which(all.ps > cuts[i-1] & all.ps <= cuts[i]);
    strata.trim[inx] <- i-1;
  }
  strata <- strata.trim
  strata[which(all.ps<cuts[1])] <- 1
  strata[which(all.ps>cuts[nstrata + 1])] <- nstrata

  # PS stratification (trim CH by range of PS in {CD,E} by stratum)
  data$strata <- strata.trim
  data <- data %>% dplyr::filter(!is.na(strata))
  temp <- data %>% group_by(strata,Z) %>% dplyr::summarise(prop=mean(Y), n=n()) %>% dplyr::mutate(var=prop*(1-prop)/n)
  point.est=mean(temp$prop[temp$Z==1] - temp$prop[temp$Z==0])
  sd = sqrt(sum(temp$var)/nstrata^2)
  delta=data.frame(point.est, lower.CI=point.est + qnorm(0.025)*sd, upper.CI=point.est + qnorm(0.975)*sd)
  outcome15 <- cbind(delta,label="PS stratification (trim)")

  comp.freq <- rbind(outcome1, outcome2, outcome3, outcome4, outcome6, outcome12, outcome15) %>%
    dplyr::mutate(sim=seed, set.num = set.num)
  return(list(mydata=mydata,comp.freq=comp.freq))
}
