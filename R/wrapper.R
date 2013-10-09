# R wrapper function for the main MCMC loop

svsample <- function(y, draws = 10000, burnin = 1000, priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1, thinpara = 1, thinlatent = 1, thintime = 1, quiet = FALSE, startpara, startlatent, expert, ...) {
 
 # Some error checking for y
 if (is(y, "svsim")) {
  y <- y[["y"]]
  warning("Extracted data vector from 'svsim'-object.")
 }
 if (!is.numeric(y)) stop("Argument 'y' (data vector) must be numeric.")
 
 if (length(y) < 2) stop("Argument 'y' (data vector) must contain at least two elements.")

 if (any(y == 0)) {
  myoffset <- sd(y)/10000
  warning(paste("Argument 'y' (data vector) contains zeros. I am adding an offset constant of size ", myoffset, " to do the auxiliary mixture sampling. If you want to avoid this, you might consider de-meaning the returns before calling this function.", sep=""))
 } else myoffset <- 0

 # Some error checking for draws
 if (!is.numeric(draws) | draws < 1) {
  stop("Argument 'draws' (number of MCMC iterations after burn-in) must be a single number >= 1.")
 } else {
  draws <- as.integer(draws)
 }

 # Some error checking for burnin
 if (!is.numeric(burnin) | burnin < 0) {
  stop("Argument 'burnin' (burn-in period) must be a single number >= 0.")
 } else {
  burnin <- as.integer(burnin)
 }
 
 # Some error checking for the prior parameters 
 if (!is.numeric(priormu) | length(priormu) != 2) {
  stop("Argument 'priormu' (mean and sd for the Gaussian prior for mu) must be numeric and of length 2.")
 }
 
 if (!is.numeric(priorphi) | length(priorphi) != 2) {
  stop("Argument 'priorphi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
 }
 
 if (!is.numeric(priorsigma) | length(priorsigma) != 1 | priorsigma <= 0) {
  stop("Argument 'priorsigma' (scaling of the chi-squared(df = 1) prior for sigma^2) must be a single number > 0.")
 }
 
 # Some error checking for thinpara
 if (!is.numeric(thinpara) | thinpara < 1) {
  stop("Argument 'thinpara' (thinning parameter for mu, phi, and sigma) must be a single number >= 1.")
 } else {
  thinpara <- as.integer(thinpara)
 }

 # Some error checking for thinlatent
 if (!is.numeric(thinlatent) | thinlatent < 1) {
  stop("Argument 'thinlatent' (thinning parameter for the latent log-volatilities) must be a single number >= 1.")
 } else {
  thinlatent <- as.integer(thinlatent)
 }

 # Some error checking for thintime
 if (!is.numeric(thintime) | thintime < 1) {
  stop("Argument 'thintime' (thinning parameter for time) must be a single number >= 1.")
 } else {
  thintime <- as.integer(thintime)
 }

 # Some error checking for expert
 if (missing(expert)) {
  para <- 3L ; parameterization <- 'GIS_C'
  mhcontrol <- -1
  gammaprior <- TRUE
  truncnormal <- FALSE
  mhsteps <- 2L
  B011 <- 10^8
  B022 <- 10^12
 } else {
  expertnames <- names(expert)
  if (!is.list(expert) | is.null(expertnames) | any(expertnames == ""))
   stop("Argument 'expert' must be a named list with nonempty names.")
  if (length(unique(expertnames)) != length(expertnames))
   stop("No duplicate elements allowed in argument 'expert'.")
  allowednames <- c("parameterization", "mhcontrol", "gammaprior", "truncnormal", "mhsteps", "proposalvar4sigmaphi", "proposalvar4sigmatheta")
  exist <- pmatch(expertnames, allowednames)
  if (any(is.na(exist)))
   stop(paste("Illegal element '", paste(expertnames[is.na(exist)], collapse="' and '"), "' in argument 'expert'.", sep=''))
  
  expertenv <- list2env(expert) 
  
  if (exists("parameterization", expertenv)) {
   parameterization <- expert[["parameterization"]]
   if (!is.character(parameterization) | is.na(parameterization)) {
    stop("Argument 'parameterization' must be either 'centered', 'noncentered', 'GIS_C', or 'GIS_NC'.")
   }
   switch(parameterization,
    centered = para <- 1L,
    noncentered = para <- 2L,
    GIS_C = para <- 3L,
    GIS_NC = para <- 4L,
    stop("Unknown parameterization. Currently you can only use 'centered', 'noncentered', 'GIS_C', and 'GIS_NC'.")
   )
  } else {
   para <- 3L ; parameterization <- 'GIS_C'
  }

  # Remark: mhcontrol < 0 means independence proposal,
  #         mhcontrol > 0 controls stepsize of log-random-walk proposal
  if (exists("mhcontrol", expertenv)) {
   mhcontrol <- expert[["mhcontrol"]]
   if (!is.numeric(mhcontrol) | length(mhcontrol) != 1)
    stop("Argument 'mhcontrol' must be a single number.")
  } else {
   mhcontrol <- -1
  }

  # use a Gamma prior for sigma^2 in C?
  if (exists("gammaprior", expertenv)) {
   gammaprior <- expert[["gammaprior"]]
   if (!is.logical(gammaprior)) stop("Argument 'gammaprior' must be TRUE or FALSE.")
  } else {
   gammaprior <- TRUE
  }

  # use a truncated normal as proposal? (or normal with rejection step)
  if (exists("truncnormal", expertenv)) {
   truncnormal <- expert[["truncnormal"]]
   if (!is.logical(truncnormal)) stop("Argument 'truncnormal' must be TRUE or FALSE.")
  } else {
   truncnormal <- FALSE
  }
 
  if (exists("mhsteps", expertenv)) {
   mhsteps <- as.integer(expert[["mhsteps"]])
   if (mhsteps != 2L & mhsteps != 1L & mhsteps != 3L) stop("mhsteps must be 1, 2, or 3")
   if (mhsteps != 2L & mhcontrol >= 0)
    stop("Log normal random walk proposal currently only implemented for mhsteps==2.")
   if (mhsteps != 2L & !isTRUE(gammaprior))
    stop("Inverse Gamma prior currently only implemented for mhsteps==2.")
  } else {
   mhsteps <- 2L
  }

  # prior for ridge _proposal_ (variance of sigma*phi)
  if (exists("proposalvar4sigmaphi", expertenv)) {
   B011 <- expert[["proposalvar4sigmaphi"]]
   if (!is.numeric(B011) | length(B011) != 1 | B011 <= 0)
    stop("Argument 'proposalvar4sigmaphi' must be a positive number.")
  } else {
   B011 <- 10^8
  }

  # prior for ridge _proposal_ (variance of sigma*theta)
  if (exists("proposalvar4sigmatheta", expertenv)) {
   B022 <- expert[["proposalvar4sigmatheta"]]
   if (!is.numeric(B022) | length(B022) != 1 | B022 <= 0)
    stop("Argument 'proposalvar4sigmatheta' must be a positive number.")
  } else {
   B022 <- 10^12
  }
 }

 # Some input checking for startpara
 if (missing(startpara)) {
  startpara <- list(mu = -10, phi = .9, sigma = .3) 
 } else {
  if (!is.list(startpara))
   stop("Argument 'startpara' must be a list with three elements named 'mu', 'phi', 'sigma'.")
  
  if (!is.numeric(startpara[["mu"]]))
   stop('Argument \'startpara[["mu"]]\' must exist and be numeric.')
  
  if (!is.numeric(startpara[["phi"]]))
   stop('Argument \'startpara[["phi"]]\' must exist and be numeric.')
  
  if (abs(startpara[["phi"]]) >= 1)
   stop('Argument \'startpara[["phi"]]\' must be between -1 and 1.')

  if (!is.numeric(startpara[["sigma"]]))
   stop('Argument \'startpara[["sigma"]]\' must exist and be numeric.')
  
  if (startpara[["sigma"]] <= 0)
   stop('Argument \'startpara[["sigma"]]\' must be positive.')
 }

 # Some input checking for startlatent
 if (missing(startlatent)) {
  startlatent <- rep(-10, length(y))
 } else {
  if (!is.numeric(startlatent) | length(startlatent) != length(y))
   stop("Argument 'startlatent' must be numeric and of the same length as the data 'y'.")
 }

 if (!quiet) {
  cat(paste("\nCalling ", parameterization, " MCMC sampler with ", draws+burnin, " iter. Series length T = ", length(y), ".\n",sep=""), file=stderr())
  flush.console()
 }

 if (.Platform$OS.type != "unix") myquiet <- TRUE else myquiet <- quiet  # Hack to prevent console flushing problems with Windows

  runtime <- system.time(res <-
  .Call("sampler", y, draws, burnin,
       	priormu[1], priormu[2]^2, priorphi[1], priorphi[2], priorsigma, 
       	thinlatent, thintime, startpara, startlatent, myquiet, para,
	mhsteps, B011, B022, mhcontrol, gammaprior, truncnormal,
	myoffset, PACKAGE = "stochvol"))

 if (any(is.na(res))) stop("Sampler returned NA. This is most likely due to bad input checks and shouldn't happen. Please report to package maintainer.")
  
 if (!quiet) {
  cat("Timing (in seconds):\n", file=stderr())
  print(runtime)
  cat(round((draws+burnin)/runtime[3]), "iterations per second.\n\n", file=stderr())
  cat("Converting results to coda objects... ", file=stderr())
 }
  
 # store results:
 # remark: +1, because C-sampler also returns the first value
 res$y <- y
 res$para <- mcmc(res$para[seq(burnin+thinpara+1, burnin+draws+1, thinpara),,drop=FALSE], burnin+thinpara, burnin+draws, thinpara)
 res$latent <- mcmc(t(res$latent), burnin+thinlatent, burnin+draws, thinlatent)
 res$latent0 <- mcmc(res$latent0, burnin+thinlatent, burnin+draws, thinlatent)
 attr(res$para, "dimnames") <- list(NULL, c("mu", "phi", "sigma"))
 attr(res$latent, "dimnames") <- list(NULL, paste('h_', seq(1, length(y), by=thintime), sep=''))
 res$runtime <- runtime
 res$priors <- list(mu = priormu, phi = priorphi, sigma = priorsigma)
 res$thinning <- list(para = thinpara, latent = thinlatent, time = thintime)
 class(res) <- "svdraws"
 
 if (!quiet) {
  cat("Done!\n", file=stderr())
  cat("Summarizing posterior draws... ", file=stderr())
 }
 res <- updatesummary(res, ...)

 if (!quiet) cat("Done!\n\n", file=stderr())
 res
}

# This function does not check input nor converts the result to coda objects

.svsample <- function(y, draws = 1, burnin = 0, priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1, thinpara = 1, thinlatent = 1, thintime = 1, quiet = TRUE, startpara, startlatent) {

 res <- .Call("sampler", y, draws, burnin, priormu[1], priormu[2]^2,
	      priorphi[1], priorphi[2], priorsigma, thinlatent, thintime,
	      startpara, startlatent, quiet, 3L, 2L, 10^8, 10^12,
	      -1, TRUE, FALSE, 0, PACKAGE = "stochvol")

 res$para <- t(res$para[-1,,drop=FALSE])
 rownames(res$para) <- names(res$para) <- c("mu", "phi", "sigma")
 res
}
