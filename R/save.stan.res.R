#' Save simulation results from Bayesian-based methods (PS + Bayesian)
#'
#' @param simu.server.path path of simulation results on server.
#' @param filename name of file for summarized simulation results, ended with ".Rdata". For example, "stan_5000.RData".
#' @param res.path path to save summarized simulation results, ended with "\\" or "/".
#' @param nsim simulation number.
#'
save.stan.res <- function(simu.server.path,
                          filename,
                          res.path,
                          nsim){
  setwd(simu.server.path)
  files <- list.files(path = simu.server.path, pattern = ".Rdata")
  comps <- list()
  print(length(files))
  num.files <- min(nsim, length(files))
  for(iter in 1:num.files){
    load(files[iter]) #comps.all
    comps[[iter]] <- comps.all
  }
  comps <- rbindlist(comps)

  save(comps, file = paste0(res.path, filename))
}
