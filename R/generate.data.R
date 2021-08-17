#' Generate data for simulation
#'
#' This function allows you to generate data for simulation
#' @param seed a single value, interpreted as an integer, or NULL. seed of random number generator. Default to 2020
#' @param n.E a positive integer for number of patients in current experimental arm (E). Default to 80.
#' @param n.CD a positive integer for number of patients in current control arm (CD). Default to 40.
#' @param n.CH a positive integer for number of patients in historical control arm (CH). Default to 300.
#' @param lambda log odds of treatment effect. HR=exp(lambda).
#' @param xi log odds of current study effect. = 0 if no drift.
#' @param x1.mean.ECD mean of continuous variable X1 in the current trial. Default to 40
#' @param x1.mean.CH mean of continuous variable X1 in the historical trial. Default to 40
#' @param x1.sd standard deviation of continuous variable X1 in both the current and historical trials. Default to 5.
#' @param x2.p.ECD mean of binary variable X2 in the current trial. Default to 0.4
#' @param x2.p.CH mean of binary variable X2 in the historical trial. Default to 0.4
#' @param x3.sd standard deviation of continuous variable X3 in both the current and historical trials. Default to 1.
#' @param beta0 intercept in the logistic outcome model. Default to 3.6
#' @param beta1 coefficient of X1 in the logistic outcome model. Default to -0.1
#' @param beta2 coefficient of X2 in the logistic outcome model. Default to -0.5
#' @param beta3 coefficient of X3 in the logistic outcome model. Default to 1
#' @param beta4 coefficient of X4 in the logistic outcome model. Default to -1
#' @keywords AugRCT
#' @import dplyr
#' @details Returns mydata, a data frame of simulated dataset, including ID, Z (indicator of treatment), D (indicator of current trial), group ("E", "CD", "CH"), x1, x2, x3, x4, Y (binary outcome)
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu
#' @export


generate.data <- function(n.E=80,n.CD=40,n.CH=300,
                          x1.mean.ECD=40,x1.mean.CH=40,x1.sd=5,
                          x2.p.ECD=0.4,x2.p.CH=0.4,
                          x3.sd=1,
                          lambda=1.19,
                          beta0=3.6,beta1=-0.1,beta2=-0.5,beta3=1,beta4=-1,
                          xi=0,
                          seed=2020){
  set.seed(seed)
  #Step 1: generate covariates and treatment indicators.
  Z <- c(rep(1, n.E), rep(0, n.CD + n.CH))
  ID <- 1: (n.E + n.CD + n.CH)
  H <- c(rep(0,n.E + n.CD),rep(1, n.CH)) #indicator of historical trial. manimulate it to allow henerogeineity in data source for CH
  D <- 1-H #indicator of current trial
  group <- c(rep("E",n.E),rep("CD",n.CD),rep("CH",n.CH))
  x1 <- c(rnorm(n = n.E + n.CD, mean = x1.mean.ECD, sd= x1.sd),rnorm(n.CH, mean = x1.mean.CH, sd = x1.sd)) # age: continuous
  x2 <- c(rbinom(n = n.E + n.CD,size = 1, prob = x2.p.ECD),rbinom(n = n.CH, size = 1,prob = x2.p.CH))# ECOG (Eastern Cooperative Oncology Group) (0-1; binary)
  x3 <- rnorm(n = n.E + n.CD + n.CH, mean = 0, sd = x3.sd) # Lab variable (continuous)
  x4 <- rbinom(n = n.E + n.CD + n.CH,size = 1, prob = 0.4)# biomarker (binary)

  #step 2: generate outcome. logistic regression. Y ~ X1 + X2 + X3 + X4
  Y <- rbinom(n = n.E + n.CD + n.CH,size = 1,
              prob = 1/(1+exp(-(beta0 + beta1*x1 + beta2*x2 + beta3* x3 + beta4*x4 + lambda*Z + xi*H))))# indicator of response

  mydata <- data.frame(ID, Z, D, group, x1, x2, x3, x4, Y) %>% dplyr::mutate(CD.ind=1*(group=="CD"),CH.ind=1*(group=="CH"))
  return(mydata)
}
