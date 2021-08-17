#' Extract the parameter of interest from a list of fits of Bayesian models (hybrid approaches)
#'
#' This function allows you to exract posterior estimates of interest from a list of fits of Bayesian models (hybrid approaches)
#' @param fits a vector of rstan fits
#' @param labels a vector of character labels for fits of hybrid approaches
#' @param par the parameter of interest you want to extract from the fits. Default to "delta".
#' @return res is the dataset with extracted estimates from fits
#' @keywords AugRCT
#' @details Returns res, a data frame of the point estimate and percentiles of a parameter of interest ( default is response rate difference between experiment and
#' augmented control group) of all the hybrid methods.
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu

Extract.Hybrid <- function(fits, labels, par="delta"){
  res <- NULL
  for (i in 1:length(labels)){
    label = labels[i]
    fit = fits[i]
    temp <- data.frame(summary(fit[[1]], pars=par)$summary, fit=label)
    res <- rbind(res, temp)
  }
  return(res)
}
