#' Save simulation results from frequentist methods (PS only methods)
#'
#' @param simu.server.path path of simulation results on server.
#' @param filename name of file for summarized simulation results, ended with ".Rdata". For example, "freq_5000.RData".
#' @param res.path path to save summarized simulation results, ended with "\\" or "/".
#' @param nsim simulation number.
#'
save.freq.res <- function(simu.server.path,
                          filename,
                          res.path,
                          nsim){
  setwd(simu.server.path)
  files <- list.files(path = simu.server.path, pattern = ".Rdata")
  comps <- list()
  print(paste0("Number of real simulations:", length(files)))
  num.files <- min(nsim, length(files)) # there might be convergence failure in some simulation
  for(iter in 1:num.files){
    load(files[iter]) #comps.all
    comps[[iter]] <- comps.freq.all
  }
  comps.freq <- rbindlist(comps)

  save(comps.freq, file = paste0(res.path, filename))
}
