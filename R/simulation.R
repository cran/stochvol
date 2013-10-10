svsim <- function(len, mu = -10, phi = 0.98, sigma = 0.2) {
 
 # Some error checking
 if (any(is.na(len)) | !is.numeric(len) | length(len) != 1 | any(len < 1)) {
  stop("Argument 'len' (length of simulated series) must be a single number >= 2.")
 } else {
  len <- as.integer(len)
 }

 if (!is.numeric(mu) | length(mu) != 1) {
  stop("Argument 'mu' (level of latent variable) must be a single number.")
 }

if (!is.numeric(phi) | length(phi) != 1) {
  stop("Argument 'phi' (persistence of latent variable) must be a single number.")
 }

if (!is.numeric(sigma) | length(sigma) != 1 | sigma <= 0) {
  stop("Argument 'sigma' (volatility of latent variable) must be a single number > 0.")
 }

 h <- rep(as.numeric(NA), len)
 h0 <- rnorm(1, mean=mu, sd=sigma/sqrt(1-phi^2))
 nu <- rnorm(len)

 # simulate w/ simple loop
 h[1] <- mu + phi*(h0-mu) + sigma*nu[1]
 for (i in seq(2, len = len-1)) h[i] <- mu + phi*(h[i-1]-mu) + sigma*nu[i]

 y <- rnorm(len, mean=0, sd=sqrt(exp(h)))  # "log-returns"
 ret <- list(y = y, vol = exp(h/2), vol0 = exp(h0/2),
	     para = list(mu = mu,
			 phi = phi,
			 sigma = sigma))
 class(ret) <- "svsim"
 ret
}

print.svsim <- function(x, ...) {
 cat("\nSimulated time series consisting of", length(x$y), "observations.\n
Parameters: level of latent variable                  mu =", x$para$mu, "
	    persistence of latent variable           phi =", x$para$phi, "
	    standard deviation of latent variable  sigma =", x$para$sigma, "
	    ")
 cat("\nSimulated initial volatility:", x$vol0, "\n")
 cat("\nSimulated volatilities:\n")
 print(x$vol, ...)
 cat("\nSimulated data (usually interpreted as 'log-returns'):\n")
 print(x$y, ...)
}

plot.svsim <- function(x, mar = c(3, 2, 2, 1), mgp = c(1.8, .6, 0), ...) {
 op <- par(mfrow = c(2, 1), mar = mar, mgp = mgp)
 plot.ts(100*x$y, ylab = "", ...)
 mtext("Simulated data: 'log-returns' (in %)", cex = 1.2, line = .4, font = 2)
 plot.ts(100*x$vol, ylab = "", ...)
 mtext("Simulated volatilities (in %)", cex = 1.2, line = .4, font = 2)
 par(op)
}

summary.svsim <- function(object, ...) {
 ret <- vector("list")
 class(ret) <- "summary.svsim"
 ret$len <- length(object$y)
 ret$para <- object$para
 ret$vol0 <- 100*object$vol0
 ret$vol <- summary(100*object$vol)
 ret$y <- summary(100*object$y)
 ret
}

print.summary.svsim  <- function(x, ...) {
 cat("\nSimulated time series consisting of ", x$len, " observations.\n",
     "\nParameters: level of latent variable                  mu = ",
     x$para$mu, 
     "\n            persistence of latent variable           phi = ",
     x$para$phi,
     "\n            standard deviation of latent variable  sigma = ",
     x$para$sigma, "\n", sep="")

 cat("\nSimulated initial volatility (in %): ")
 cat(x$vol0, "\n")
 cat("\nSummary of simulated volatilities (in %):\n")
 print(x$vol)
 cat("\nSummary of simulated data (in %):\n")
 print(x$y)
 invisible(x)
}
