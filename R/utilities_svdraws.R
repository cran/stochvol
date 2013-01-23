para <- function(x) {
 x$para
}

latent <- function(x) {
 x$latent
}

latent0 <- function(x) {
 x$latent0
}

priors <- function(x) {
 x$priors
}

thinning <- function(x) {
 x$thinning
}

runtime <- function(x) {
 x$runtime
}

updatesummary <- function(x, quantiles = c(.05, .5, .95), esspara = TRUE, esslatent = FALSE) {
 
 summaryfunction <- function(x, quants = quantiles, ess = TRUE) {
  if (ess) {
   c(mean = mean(x), sd = sd(x), quantile(x, quantiles),
    ESS = as.numeric(effectiveSize(x)))
  } else {
   c(mean = mean(x), sd = sd(x), quantile(x, quantiles))
  }
 }

 res <- list()
 
 res$para <- t(apply(x$para, 2, summaryfunction, ess = esspara))
 res$para <- rbind(res$para, "exp(mu/2)" = c(summaryfunction(exp(x$para[,"mu"]/2), ess=FALSE), res$para["mu", "ESS"]))
 res$para <- rbind(res$para, "sigma^2" = c(summaryfunction(x$para[,"sigma"]^2, ess=FALSE), res$para["sigma", "ESS"]))
 
 res$latent <- t(apply(x$latent, 2, summaryfunction, ess = esslatent))
 tmp <- exp(x$latent/2)
 res$latent <- cbind(res$latent, "mean(exp(h_t/2))" = colMeans(tmp))
 res$latent <- cbind(res$latent, "sd(exp(h_t/2))" = apply(tmp, 2, sd))
 
 res$latent0 <- c(summaryfunction(x$latent0, ess = esslatent), "mean(exp(h_t/2))" = mean(exp(x$latent0/2)), "sd(exp(h_t/2))" = sd(exp(x$latent0/2)))
 
 x$summary <- res
 x
}

summary.svdraws <- function(object, showpara = TRUE, showlatent = TRUE, ...) {
 mcp <- mcpar(para(object))
 mcl <- mcpar(latent(object))

 cat("\nSummary of ", mcp[2]-mcp[1]+mcp[3], " MCMC draws after a burn-in of ", mcp[1]-mcp[3], ".\n
Prior distributions:
mu        ~ Normal(mean = ", priors(object)$mu[1], ", sd = ", priors(object)$mu[2], ")
(phi+1)/2 ~ Beta(a0 = ", priors(object)$phi[1], ", b0 = ", priors(object)$phi[2], ")
sigma^2   ~ ", priors(object)$sigma, " * Chisq(df = 1)\n", sep="")

 if (all(isTRUE(showpara))) {
  cat("\nPosterior draws of parameters (thinning = ", mcp[3], "):\n", sep='')
  print(para(object$summary), digits=2, ...)
 }
 
 if (all(isTRUE(showlatent))) {
  cat("\nPosterior draws of initial and contemporaneous latents (thinning = ", mcl[3], "):\n", sep='')
  print(rbind("h_0" = latent0(object$summary), latent(object$summary)), digits=2, ...)
 }
 invisible(object)
}

print.svdraws <- function(x, showpara = TRUE, showlatent = TRUE, ...) {
 if (all(isTRUE(showpara))) {
  cat("\n*** Posterior draws of parameters ***\n")
  print(para(x), ...)
 }
 
 if (all(isTRUE(showlatent))) {
  cat("\n*** Posterior draws of initial latent variable h_0 ***\n")
  print(latent0(x), ...)
  cat("\n*** Posterior draws of contemporaneous latent variables h_t ***\n")
  print(latent(x), ...)
 }
 invisible(x)
}

paradensplot <- function(x, showobs = TRUE, showprior = TRUE, showxlab = TRUE,
		       	 mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), ...) {
 if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
 oldpar <- par(mar=mar)
 paranames <- c(quote(mu), quote(phi), quote(sigma))
 cutat1 <- c(FALSE, TRUE, FALSE)
 for (i in 1:3) {
  mydensplot(x$para[,i], show.obs=showobs, main=paste("Density of", paranames[i]),
	     cutat1=cutat1[i], showxlab=showxlab, mgp = mgp, ...)
  if (isTRUE(showprior)) {
   paras <- x$priors[[i]]
   vals <- seq(from=par('usr')[1], to=par('usr')[2], len=1000)
   if (i == 1) lines(vals, dnorm(vals, paras[1], paras[2]), col=8, lty=2)
   else if (i == 2) lines(vals, .5*dbeta((vals+1)/2, paras[1], paras[2]), col=8, lty=2)
   else if (i == 3) lines(vals, 2*dnorm(vals, 0, sqrt(paras[1])), col=8, lty=2)
  }
 }
 par(oldpar)
 invisible(x)
}

paratraceplot <- function(x, mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), ...) {
 if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
 oldpar <- par(mar=mar)
 paranames <- c(quote(mu), quote(phi), quote(sigma))
 for (i in 1:3) {
  mytraceplot(x$para[,i], xlab="", mgp = c(2, .6, 0),
	      main=paste0("Trace of ", paranames[i], " (thinning = ", x$thinning$para,")"), ...)
 }
 par(oldpar)
 invisible(x)
}

volplot <- function(x, forecast = 0, dates = NULL, show0 = FALSE,
		    col = NULL, forecastlty = NULL, tcl = -.4,
		    mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0), ...) {
 if (!is(x, "svdraws")) stop("This function expects an 'svdraws' object.")
 oldpar <- par(mar = mar)
 where <- grep("%", dimnames(x$summary$latent)[[2]])
 obj <- t(100*exp(x$summary$latent[,where]/2))  # monotone transformation!
 qs <- dim(obj)[1]
 timelen <- dim(obj)[2]
 if (is.null(qs) | all(is.na(obj))) stop("No quantiles to plot.")
 if (is.null(col)) {
  cols <- rep(8, qs)
  cols[dimnames(obj)[[1]] == "50%"] <- 1
 } else cols <- col
 if (is.null(forecastlty)) forecastlty <- 2
 
 if (is(forecast, "svpredict") | (is.numeric(forecast) & length(forecast) == 1 & all(forecast != 0))) { # also draw future values
  thintime <- x$thinning$time
  
  if (thintime != 1) {
   lasth <- as.integer(gsub("h_", "", dimnames(x$latent)[[2]][dim(x$latent)[2]]))
   if (length(x$y) > lasth) {
       warning(paste0("Thinning for time 'thintime' has not been set to one during sampling. This means we are forecasting conditional on h_", lasth, " and not on h_", length(x$y), "."))
   }
  }
  
  if(is.numeric(forecast) & length(forecast) == 1 & all(forecast >= 1)) {
#   warning("Calling prediction method.")
   forecast <- predict(x, forecast)
  }
  if(!is(forecast, "svpredict")) stop("Argument 'forecast' must be a single nonnegative integer, or of class type 'svpredict'.")
  
  futlen <- dim(forecast)[2]
  
  xs <- matrix(rep(seq(timelen, timelen + futlen/thintime, len=futlen+1), qs), nrow=futlen+1)
  quants <- as.numeric(gsub("%", "", dimnames(obj)[[1]]))/100
  ys <- rbind(obj[,timelen], t(matrix(apply(100*exp(forecast/2), 2, quantile, quants), nrow=qs)))
 
  if (futlen/thintime > .01*timelen) {  # increase xlim to give space for forecast
   if (thintime == 1) {
    xlim <- c(0, timelen+futlen/thintime)
   } else {
    xlim <- c(1, timelen+futlen/thintime)
   }
  } else {
   xlim <- NULL
  }
 } else xlim <- NULL
 
 ts.plot(t(obj), gpars=list(xlim=xlim, col=cols, xlab='', xaxt='n', mgp=mgp, tcl=tcl,
			    main = paste0("Estimated volatilities in percent (",
					  paste(dimnames(obj)[[1]], collapse=' / '),
					  " posterior quantiles)"), ...))
 
 if (is(forecast, "svpredict")) {
  for (i in 1:qs) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
 }
 
 ax <- axis(1, tick=FALSE, labels=FALSE)  # just automagic axis ticks, don't draw yet

 if (show0) { # also draw latent0:
  thintime <- x$thin$time
  xs <- matrix(rep(c(1-1/thintime,1), qs), nrow=2)
  where <- grep("%", names(x$summary$latent0))
  ys <- rbind(100*exp(x$summary$latent0[where]/2), obj[,1])
  for (i in 1:qs) lines(xs[,i], ys[,i], lty=forecastlty, col=cols[i])
 }
   
     
 if (is.null(dates)) {
  dates <- c(0L, as.integer(gsub("h_", "", dimnames(x$latent)[[2]])))
  if (max(ax) > length(dates)) {  # means we are probably forecasting and need extra axis labels
   dates <- c(dates, seq(length(dates), max(ax), by=dates[2]-dates[1]))
  }
 }
 else {
  if (length(dates) != ncol(x$latent)) {
   warning("Length of argument 'dates' differs from ncol(x$latent).")
  }
  dates <- c('', dates)
  ax <- ax[ax != 0]  # avoid "zero" tick
 }
 axis(1, at=ax, labels=dates[ax+1], mgp=mgp, tcl=tcl)
 par(oldpar)
 invisible(x)
}

plot.svdraws <- function(x, forecast = NULL, dates = NULL, show0 = FALSE,
			 showobs = TRUE, showprior = TRUE, col = NULL,
			 forecastlty = NULL, tcl = -0.4,
			 mar = c(1.9, 1.9, 1.7, .5), mgp = c(2, .6, 0),
			 ...) {
 layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), 3, 3, byrow = TRUE))
 volplot(x, dates = dates, show0 = show0, forecast = forecast,
	 forecastlty = forecastlty, col = col, tcl = tcl, mar = mar,
	 mgp = mgp, ...)
 paratraceplot(x, mar = mar, mgp = mgp, ...)
 paradensplot(x, showobs = showobs, showprior = showprior,
	      showxlab = FALSE, mar = mar, mgp = mgp, ...)
 invisible(x)
}

predict.svdraws <- function(object, steps = 1, ...) {
 if (!is(object, "svdraws")) stop("Argument 'object' must be of class 'svdraws'.")
 steps <- as.integer(steps)
 if (steps < 1) stop("Argument 'steps' must be greater or equal to 1.")
 thinlatent <- object$thinning$latent
 thinpara <- object$thinning$para
 if (thinpara != thinlatent) {
  warning("Thinning of parameters is different from thinning of latent variables. Trying to sort this out.")
  if (thinpara %% thinlatent == 0) {
   usepara <- 1:(dim(object$para)[1])
   uselatent <- seq(thinpara/thinlatent, dim(object$latent)[1], by=thinpara/thinlatent)
  } else if (thinlatent %% thinpara == 0) {
   uselatent <- 1:(dim(object$latent)[1])
   usepara <- seq(thinlatent/thinpara, dim(object$para)[1], by=thinlatent/thinpara)
  } else stop("Incompatible thinning parameters. Prediction currently not implemented.")
 } else {
  usepara <- uselatent <- seq.int(dim(object$para)[1])
 }
 
 mu <- object$para[,"mu"][usepara]
 phi <- object$para[,"phi"][usepara]
 sigma <- object$para[,"sigma"][usepara]
 hlast <- object$latent[,dim(object$latent)[2]][uselatent]

 mythin <- max(thinpara, thinlatent)
 len <- length(sigma)
 volpred <- mcmc(matrix(as.numeric(NA), nrow=len, ncol=steps), start=mythin, end=len*mythin, thin=mythin)
 
 volpred[,1] <- mu + phi*(hlast - mu) + sigma*rnorm(len)
 if (steps > 1)
  for (i in (seq.int(steps-1) + 1))
   volpred[,i] <- mu + phi*(volpred[,i-1] - mu) + sigma*rnorm(len)
 
 class(volpred) <- c("svpredict", "mcmc")
 lastname <- dimnames(object$latent)[[2]][dim(object$latent)[2]]
 lastnumber <- as.integer(gsub("h_", "", lastname))
 colnames(volpred) <- paste0("h_", seq(lastnumber + 1, lastnumber + steps))
 volpred
}


# modified density plot (from coda package)
mydensplot <- function (x, show.obs = TRUE, bwf, main = "", ylim, cutat1=FALSE, showxlab=TRUE, mgp = c(2,.6,0), tcl=-.4, ...) 
{
    xx <- as.matrix(x)
    for (i in 1:nvar(x)) {
        y <- xx[, i, drop = TRUE]
        if (missing(bwf)) 
            bwf <- function(x) {
                x <- x[!is.na(as.vector(x))]
                return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
            }
        bw <- bwf(y)
        width <- 4 * bw
        if (max(abs(y - floor(y))) == 0 || bw == 0) 
            hist(y, prob = TRUE, main = main, ...)
        else {
            scale <- "open"
        if (isTRUE(cutat1)) {
	     if (1-max(y) < 2*bw) {
	      scale <- "cutat1"
	      y <- c(y, 2 - y)
	      if (1+min(y) < 2*bw) {
	       scale <- "cutatboth"
	       y <- c(y, -2 - y, 2 - y)
	      }
	     } else if (1+min(y) < 2*bw) {
	      scale <- "cutat-1"
	      y <- c(y, -2 - y)
	     }
	    }
	    else if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
            	if (min(y) >= 0 && min(y) < 2 * bw) {
                  scale <- "proportion"
                  y <- c(y, -y, 2 - y)
                }
            }
            else if (min(y) >= 0 && min(y) < 2 * bw) {
                scale <- "positive"
                y <- c(y, -y)
            }
	    else scale <- "open"
            dens <- density(y, width = width)
            if (scale == "proportion") {
                dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 
                  1]
                dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
            }
            else if (scale == "positive") {
                dens$y <- 2 * dens$y[dens$x >= 0]
                dens$x <- dens$x[dens$x >= 0]
            }
	    else if (scale == "cutat1") {
	    	dens$y <- 2 * dens$y[dens$x <= 1]
	    	dens$x <- dens$x[dens$x <= 1]
	    }
	    else if (scale == "cutat-1") {
	    	dens$y <- 2 * dens$y[dens$x >= -1]
	    	dens$x <- dens$x[dens$x >= -1]
	    }
	    else if (scale == "cutatboth") {
	    	dens$y <- 3 * dens$y[dens$x >= -1 & dens$x <= 1]
	    	dens$x <- dens$x[dens$x >= -1 & dens$x <= 1]
	    }
	    if (missing(ylim)) 
                ylim <- c(0, max(dens$y))
            
	    plot(dens, ylab = "", main = main, type = "l", 
		  ylim = ylim, xlab="", mgp = mgp, tcl = tcl, ...)
            if(isTRUE(showxlab)) {
	       if (is.R()) {
                  mtext(paste("N =", niter(x), "  Bandwidth =",
			      formatC(dens$bw)), side=1, line=2.7, cex=.7)
               } else {
                  mtext(paste("N =", niter(x), "  Bandwidth =",
			      formatC(bw)), side=1, line=2.7, cex=.7)
               }
	    }
            if (show.obs) 
                lines(y[1:niter(x)], rep(max(dens$y)/100, niter(x)), 
                  type = "h")
        }
        if (!is.null(varnames(x)) && is.null(list(...)$main)) 
            title(paste("Density of", varnames(x)[i]))
    }
    return(invisible(x))
}

# modified traceplot (from coda)
mytraceplot <- function (x, smooth = FALSE, col = 1:6, type = "l", ylab = "", xlab = "Iterations", mgp = c(2,.6,0), tcl = -.4, ...) 
{
    x <- mcmc.list(x)
    args <- list(...)
    for (j in 1:nvar(x)) {
        xp <- as.vector(time(x))
        yp <- if (nvar(x) > 1) 
            x[, j, drop = TRUE]
        else x
        yp <- do.call("cbind", yp)
        matplot(xp, yp, xlab = xlab, ylab = ylab, type = type, 
            col = col, mgp = mgp, tcl = tcl, ...)
        if (!is.null(varnames(x)) && is.null(list(...)$main)) 
            title(paste0("Trace of ", varnames(x)[j], " (thin = ", attr(x, "thinning")$thinpara,")"))
        if (smooth) {
            scol <- rep(col, length = nchain(x))
            for (k in 1:nchain(x)) lines(lowess(xp, yp[, k]), 
                col = scol[k])
        }
    }
}
