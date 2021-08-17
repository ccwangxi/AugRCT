#' Implement Frequentist (naive and PS only) approaches in real data application
#'
#' This function allows you to implement all the frequentist approaches (naive and PS only) in augmenting current control with historical data in real data application.
#' @param realdata the name of a data frame containing all needed variables, including D (indicator of current trial), Z (indicator of treatment),
#' Y (binary outcome), a list of baseline covariates included in PS model
#' @param ps.covs a list of characters for variable na,es included in PS model. e.g. ps.covs=c("X1","X2")
#' @param trt a character of variables name for the binary outcome in PS model. Defaul to "D"
#' @param nstrata a positive integer of number of strata in PS stratification. Default to 5.
#' @param seed a single value, interpreted as an integer, or NULL. seed of random number generator. Default to 2020
#' @param bootstrap.seed seed for bootstrap CI of IPTW approach. Default to 1234
#' @param bootstrap.rep replication number of Bootstrap samples. Default to 10000
#' @return comp.freq, a data frame of the point estimate and percentiles of delta (response rate difference between experiment and augmented control group) of all the frequentist methods.
#' @return pm.realdata, a data frame with PS matching adjusted patients.
#' @return realdata.iptw, a data frame with IPTW adjusted patients.
#' @return strat.trim.realdata, a data frame with PS stratification adjusted patients. variable strata in the data frame
#' indicates stratum number of each subject.
#' @keywords AugRCT
#' @import dplyr tidyr MatchIt rstan
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu
#' @export

RealData.Freq <- function(realdata,
                          ps.covs, nstrata = 5,
                          trt="D", seed=2020,
                          bootstrap.seed=1234, bootstrap.rep=10000){
  # this funciton is used to implement the available methods (Freq) to the real data
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
  pm.realdata.check <- realdata1 %>% #dplyr::filter(group=="CD" | (!is.na(matches))) %>%
    dplyr::mutate(keep=ifelse(group %in% c("E","CD"),1,
                              ifelse(is.na(matches), 0, NA)))
  random.list <- sample(c(rep(0,drop.num),rep(1,keep.num)),replace=FALSE)
  pm.realdata.check$keep[is.na(pm.realdata.check$keep)] <- random.list
  pm.realdata <- pm.realdata.check %>%
    dplyr::filter(keep==1)

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

  #------------------------------------------------
  # Frequentist approaches
  #PS Matching + logistic outcome model
  outcome1 <- get.mn.ci(data = pm.realdata, label = "Matching")

  # IPTW (ATT)  + logistic outcome model (MSM marginal structure model. no truncation of weights for now)
  prop <- realdata %>% dplyr::group_by(D) %>% dplyr::summarise(mean = mean(Y)); prop.D = prop$mean[prop$D==1]
  suppressMessages({realdata.iptw <- realdata %>%
    dplyr::mutate(weight.ns = ifelse(realdata$group %in% c("E","CD"),1, est.PS/(1-est.PS)))})# weight.ns: non-stabilized weight
  realdata.iptw$weight[realdata.iptw$D==0] <- realdata.iptw$weight.ns[realdata.iptw$D==0]/mean(realdata.iptw$weight.ns[realdata.iptw$D==0])
  if (max(realdata.iptw$weight) > 10) {print("Stabilized IPTW weight > 10")}

  data <- realdata.iptw
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
    data <- realdata.iptw[x,]
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
  outcome3 <- get.mn.ci(data = realdata, label = "Pool")
  # logistic outcome model (no borrowing)
  outcome4 <- get.mn.ci(data = realdata %>% dplyr::filter(group != "CH"), label = "No borrowing")

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
  data <- fit.trim.1(data=realdata)
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
    data <- realdata.iptw[x,]
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
  data <- realdata
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
  strat.trim.realdata <- data
  data <- data %>% dplyr::filter(!is.na(strata))
  temp <- data %>% group_by(strata,Z) %>% dplyr::summarise(prop=mean(Y), n=n()) %>% dplyr::mutate(var=prop*(1-prop)/n)
  point.est=mean(temp$prop[temp$Z==1] - temp$prop[temp$Z==0])
  sd = sqrt(sum(temp$var)/nstrata^2)
  delta=data.frame(point.est, lower.CI=point.est + qnorm(0.025)*sd, upper.CI=point.est + qnorm(0.975)*sd)
  outcome15 <- cbind(delta,label="PS stratification (trim)")

  comp.freq <- rbind(outcome1, outcome2, outcome3, outcome4, outcome12, outcome15)
  return(list(pm.realdata = pm.realdata,
              comp.freq = comp.freq,
              realdata.iptw = realdata.iptw,
              strat.trim.realdata = strat.trim.realdata,
              rwe.data = rwe.data))
}
