#' Calculate the point estimate and 95\% Bootstrap CI with IPTW approach
#'
#' This function allows you to calculate the point estimate and 95\% Bootstrap CI with IPTW approach of treatment effect in RD using IPTW1 estimator in Lunceford's paper (Stratification and Weighting Via the Propensity Score in Estimation of Causal Treatment Effects: A Comparative Study)
#' @param data the name of a data frame containing all needed variables, including Z (indicator of treatment) and Y (binary outcome).
#' @param alpha level of Type I error. Default to 0.05 to require 95\% CI.
#' @param label a character of labels for a hybrid approach
#' @keywords AugRCT
#' @import dplyr
#' @details Returns a data frame with one row, and four columns for point.est, lower.CI, upper.CI and label.
#' augmented control group) of all the frequentist methods.
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu

get.iptw.ci <- function(data,label, alpha=0.05){
  point.est
  temp <- data %>% group_by(Z) %>% dplyr::summarise(Ysum=sum(Y),n=n())
  ci <- c(diffscoreci(temp$Ysum[temp$Z==1],temp$n[temp$Z==1],temp$Ysum[temp$Z==0],temp$n[temp$Z==0], 1-alpha)$conf.int)
  res=data.frame(point.est=temp$Ysum[temp$Z==1]/temp$n[temp$Z==1] - temp$Ysum[temp$Z==0]/temp$n[temp$Z==0], lower.CI=ci[1], upper.CI=ci[2],label=label)
  res
}
