#-----------------------------------------------------
# Functions to find true values of each simulation setting
#-----------------------------------------------------
#' Find the true value of lambda such that desired power is achieved in the standard setting via simulation
#'
#' The value of lambda is found under the standard RCT 1:1 with 80 patients in experimental arm and 80 in control arm
#' @param power desired power. Default to 0.8
#' @param n.sim number of simulation for power calculation
#' @param search.lowbound lower bound to search the true value of lambda. Default to 1.
#' @param search.uppbound upper bound to search the true value of lambda. Default to 1.3.
#' @param seed.init seed for random number generator
#' @export
#' @details Return the value of lambda

find.lambda <- function(n.sim=1000, power=0.8, seed.init = 2020, search.lowbound = 1, search.uppbound = 1.3){

  find.lambda.fn <- function(x){
    set.seed(seed.init)
    seeds <- as.integer(rnorm(n.sim,0,1000))
    pow.n <- 0
    for(i in 1:n.sim){
      seed <- seeds[i]
      mydata <- generate.data(seed=seed,n.E=80,n.CD=40,n.CH=40,lambda=x,x1.mean.ECD=40,x1.mean.CH=40,
                              x2.p.ECD=0.4,x2.p.CH=0.4,
                              beta0=3.6,beta1=-0.1,beta2=-0.5,beta3=1,beta4=-1)
      #use data with n.E=80, n.C= n.CD+n.CH=80 (here, CH and CD are from the same distribution)
      fit <- glm(Y ~ Z, data = mydata, family = "binomial") #default link is logit link
      suppressMessages({pow.n <- pow.n + 1*(data.frame(confint(fit))["Z","X2.5.."] > 0)})
    }
    pow.n/n.sim - power
  }
  # get the value of lambda
  lambda <- uniroot(find.lambda.fn, tol=1e-3, interval = c(search.lowbound,search.uppbound))$root
  lambda
}

#' Find the true value of xi such that desired drift is achieved via simulation
#'
#' The value of xi is found for desired value of drift by specifying thetaCH with different treatment effect (lambda)
#' @param thetaCH desired thetaCH.
#' @param lambda log odds of treatment effect. HR=exp(lambda). Default to 1.185, the scenario with a positive treatment effect s.t. a power of 0.8 is achived in RCT 1:1 (80 patients in E group, 80 in C group)
#' @param search.lowbound lower bound to search the true value of xi. Default to -10.
#' @param search.uppbound upper bound to search the true value of xi. Default to 10.
#' @param seed seed for random number generator
#' @export
#' @details Return the value of xi

find.xi <- function(thetaCH, seed = 2020, lambda=1.185, search.lowbound = -10, search.uppbound = 10){
  # this funciton is used to find the true value of xi such that desired thetaCH is achieved
  # theta.CH: the theta.CH we want to achieve
  # lambda: for treatment effect
  # search.lowbound, search.uppbound: lower bound and upper bound to search the true value of lambda.
  find.xi.fn <- function(x){
    dt <- generate.data(seed=seed,lambda=lambda, xi=x,
                        n.E=80*10000,n.CD=40*10000,n.CH=300*10000,x1.mean.ECD=40,x1.mean.CH=40,
                        x2.p.ECD=0.4,x2.p.CH=0.4,
                        beta0=3.6,beta1=-0.1,beta2=-0.5,beta3=1,beta4=-1)
    trues <- dt %>% dplyr::group_by(group) %>% dplyr::summarise(mean=mean(Y))
    true.thetaCH <- trues$mean[trues$group=="CH"]
    true.thetaCH - thetaCH
  }
  # get the value of xi
  xi <- uniroot(find.xi.fn, tol=1e-3, interval = c(search.lowbound,search.uppbound))$root
  return(xi)
}
