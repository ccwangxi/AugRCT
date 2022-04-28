#' Manipulate results from function Simu.Hybrid
#'
#' This function allows you to compare the parameters of interest for different model fits
#' @param fits a vector of rstan fits
#' @param labels a vector of notations or names for rstan fits
#' @param par the name of parameter you want to extract
comp.fits <- function(fits, labels, par="delta"){
  comp <- NULL
  for (i in 1:length(labels)){
    label = labels[i]
    fit = fits[i]
    temp <- data.frame(summary(fit[[1]], pars=par)$summary, fit=label)
    comp <- rbind(comp, temp)
  }
  return(comp)
}
