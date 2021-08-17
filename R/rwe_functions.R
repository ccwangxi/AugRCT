########################################
# Title: R Functions for  stratification + power prior approach
# Subtitle: modify the source code from https://github.com/olssol/psrwe (Wang 2019)
# Author: Xi "Ada" Wang
# Date: June/10/2020
# Modified: July/02/2020
########################################
#' Stratify patients into nstrata
rwePS <- function(data, ps.fml = NULL,
                  v.grp = "D", #outcome variable in PS model
                  v.covs = "V1", #list of covariates
                  d1.grp = 1, #level of v.grp variable
                  nstrata = 5) {
  if (is.null(ps.fml))
    ps.fml <- as.formula(paste(v.grp, "~", paste(v.covs, collapse="+"), sep=""));

  all.ps  <- data$est.PS;
  D1.ps   <- all.ps[which(1 == data[[v.grp]])]; # CD + E

  ## stratification
  cuts    <- quantile(D1.ps, seq(0, 1,length=nstrata+1));
  cuts[1] <- cuts[1] - 0.0001;
  rst     <- rep(NA, length(all.ps));
  for (i in 2:length(cuts)) {
    inx      <- which(all.ps > cuts[i-1] & all.ps <= cuts[i]);
    rst[inx] <- i-1;
  }
  strata <- rst

  data[["_ps_"]]     <- all.ps;
  data[["_strata_"]] <- strata;
  data[["_grp_"]]    <- data$D;

  rst <- list(data    = data,
              ps.fml  = ps.fml,
              nstrata = nstrata);

  rst
}


#' Compute distance from F0 to F1
#'
#' @param type type of distances. ovl: overlapping coefficient, kl:
#'     1/(1+Kullback-Leibler divergence)
#' @param n.bins number of bins for KL computation
#' @param epsilon small integer for Dirichlet smoothing
#'
#' @return a vector with the number of samples in group 0, the number of samples
#'     in group 1, and 1/(1+KL divergence) from group 0 to group 1 when type is
#'     kl, or the overlapping coefficient when type is ovl
#'
rweDist <- function(sample.F0, sample.F1, n.bins = 10, type = c("kl", "ovl"), epsilon = 10^-6) {

  rweCut <- function(x, y=x, breaks = 5) {
    #this function is used to Cut a sequence of numbers x into bins with equal numbers in each bin
    cuts    <- quantile(x, seq(0, 1,length=breaks+1));
    cuts[1] <- cuts[1] - 0.001;
    rst     <- rep(NA, length(y));
    for (i in 2:length(cuts)) {
      inx      <- which(y > cuts[i-1] & y <= cuts[i]);
      rst[inx] <- i-1;
    }

    rst
  }

  type     <- match.arg(type);

  smps     <- c(sample.F0, sample.F1);
  n0       <- length(sample.F0);
  n1       <- length(sample.F1);

  if (1 == length(unique(smps))) {
    cut.smps <- rep(1, n0+n1)
    n.bins   <- 1;
    warning("Distributions for computing distances are degenerate.",
            call. = FALSE);
  } else {
    cut.smps <- rweCut(smps, breaks = n.bins);
  }

  rst <- 0;
  for (j in 1:n.bins) {
    n0.j <- length(which(j == cut.smps[1:n0]));
    n1.j <- length(which(j == cut.smps[(n0+1):(n0+n1)]));

    rst  <- rst + switch(type,
                         kl = {ep0  <- (n0.j+epsilon)/(n0 + epsilon * n.bins);
                         ep1  <- (n1.j+epsilon)/(n1 + epsilon * n.bins);
                         ep1 * log(ep1/ep0)},
                         ovl = min(n0.j/n0, n1.j/n1));
  }

  if ("kl" == type)
    rst <- 1/(1+rst);

  c(n0,n1,rst);
}


#' Get number of subjects and the distances of PS distributions for each PS strata
#'
#' @param data.withps data frame with estimated PS
#' @param n.bins number of bins for KL computation. Default to 10
#'
rwePSDist <- function(data.withps, n.bins = 10, type = c("ovl", "kl"), ...) {

  type <- match.arg(type);

  dataps   <- data.withps$data %>% dplyr::filter(! is.na(`_strata_`))
  nstrata  <- data.withps$nstrata;
  rst     <- NULL;
  for (i in 1:nstrata) {
    ps0 <- dataps[which(i == dataps[["_strata_"]] &
                          0 == dataps[["_grp_"]]),
                  "_ps_"];
    ps1 <- dataps[which(i == dataps[["_strata_"]] &
                          1 == dataps[["_grp_"]]),
                  "_ps_"];

    if (0 == length(ps0) | 0 == length(ps1))
      warning("No samples in strata");

    if (any(is.na(c(ps0, ps1))))
      warning("NA found in propensity scores in a strata");

    cur.dist <- rweDist(ps0, ps1, n.bins = n.bins, type = type, ...); #in each stratum, calculate overlapping coefficient, rs
    rst      <- rbind(rst, c(i, cur.dist));
  }

  ##overall
  ps0        <- dataps[which(0 == dataps[["_grp_"]]), "_ps_"];
  ps1        <- dataps[which(1 == dataps[["_grp_"]]), "_ps_"];
  all.dist   <- rweDist(ps0, ps1, n.bins = nstrata*n.bins, type = type, ...);
  rst        <- rbind(rst, c(0, all.dist)); # rs overall (i.e. all strata)


  colnames(rst) <- c("Strata", "N0", "N1", "Dist");#N0: N of trimmed CH, N1: N of CD and E
  rst           <- data.frame(rst);

  rst
}

#'  Run stan model
#'
#'  @param lst.data list of data to run STAN model
#'  @param model name of STAN model
#'  @export
rweSTAN <- function(lst.data, model,
                    chains = 4, iter = 10000, warmup = 1000,
                    control = list(adapt_delta=0.95), ...) {

  stan.rst <- rstan::sampling(model,
                              data    = lst.data,
                              chains  = chains,
                              iter    = iter,
                              warmup  = warmup,
                              control = control,
                              ...);

  stan.rst;
}

#' Get Posterior for all stratum
#'
#' @param data class DWITHPS data frame
#' @param type type of outcomes
#' @param A    target number of subjects to be borrowed
#' @param RS   parameters for dirichelet prior
#' @param ...  extra parameters for calling function \code{\link{rweSTAN}}
#'
#' @export
#'
rwePsPreprocess <- function(data, A = 0, RS = NULL,
                            v.outcome = "Y", group="group") {

  stopifnot(v.outcome %in% colnames(data));

  ## prepare stan data
  data <- data %>% dplyr::filter(! is.na(`_strata_`))
  S      <- max(data[["_strata_"]]);
  stan.d <- NULL;

  YE     <- NULL;
  INXE   <- NULL;
  for (i in 1:S) {
    cur.dE <- data[data[["_strata_"]] == i & data[,group] == "E", v.outcome]; #Y for E CD in stratum S
    cur.dCH <- data[data[["_strata_"]] == i & data[,group] == "CH", v.outcome]; #Y for CH in stratum S
    cur.dCD <- data[data[["_strata_"]] == i & data[,group] == "CD", v.outcome]; #Y for CD in stratum S

    if (0 == length(cur.dE)) {
      stop(paste("Stratum ", i, " contains no subjects from E", sep = ""));
    }

    cur.nE <- length(cur.dE);
    cur.d  <- c(NCH = length(cur.dCH), YBARCH = mean(cur.dCH), SDCH   = sd(cur.dCH), YSUMCH = sum(cur.dCH),
                NE = cur.nE, YBARE = mean(cur.dE), YSUME = sum(cur.dE),
                NCD = length(cur.dCD), YSUMCD = sum(cur.dCD));

    stan.d <- rbind(stan.d, cur.d);
    YE     <- c(YE, cur.dE);
    INXE   <- c(INXE, rep(i, length = cur.nE)); #sequence of stratum indicators for E
  }

  if (is.null(RS))
    RS <- rep(1/S, S);

  lst.data  <- list(S     = S,
                    A     = A,
                    RS    = RS,
                    NCH    = stan.d[,"NCH"], #N0 is NCH
                    NE    = stan.d[,"NE"],
                    NCD    = stan.d[,"NCD"],
                    YBARCH = stan.d[,"YBARCH"],
                    SDCH   = stan.d[,"SDCH"],
                    YSUMCH = stan.d[,"YSUMCH"]);

  ## sampling
  lst.data <- c(lst.data,
                list(YBARE = as.numeric(stan.d[,"YBARE"]),
                     YSUME = as.numeric(stan.d[,"YSUME"]),
                     YSUMCD = as.numeric(stan.d[,"YSUMCD"])));

  ## return
  lst.data;
}
