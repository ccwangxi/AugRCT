#' Add estimated propensity score to a data frame
#'
#' This function allows you to add estimated PS from propensity model logit(P("trt"=1)) = ps.covs
#' @param data the name of a data frame containing all needed variables, including D (indicator of current trial) and a list of baseline covariates included in PS model
#' @param ps.covs a list of characters for variable na,es included in PS model
#' @param trt a character of variables name for the binary outcome in PS model. Defaul to "D"
#' @keywords AugRCT
#' @details Returns data.ps, the data frame with one added column of est.PS for estimated PS.
#' @details Returns PS.fit, the glm fit of the propensity model
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu

Add.PS <- function(data,trt="D",ps.covs=c("x1","x2","x3","x4")){
  ps.fml <- as.formula(paste(trt, "~", paste(ps.covs, collapse="+"), sep=""))
  PS.fit <- glm(ps.fml, data = data, family = "binomial")
  data$est.PS <- PS.fit$fitted.values #predict(PS.logit, type="response")

  return(list(PS.fit=PS.fit,
              data.ps=data))
}
