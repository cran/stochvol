\encoding{UTF-8}
\name{stochvol-package}
\alias{stochvol-package}
\alias{stochvol}
\docType{package}
\title{
Efficient Bayesian Inference for Stochastic Volatility (SV) Models
}
\description{
This package provides an efficient algorithm for fully Bayesian estimation of stochastic volatility (SV) models via Markov chain Monte Carlo (MCMC) methods. Details can be found in Kastner and Frühwirth-Schnatter (2013).
}
\details{
Bayesian inference for stochastic volatility models using MCMC methods highly depends on actual parameter values in terms of sampling efficiency. While draws from the posterior utilizing the standard centered parameterization break down when the volatility of volatility parameter in the latent state equation is small, non-centered versions of the model show deficiencies for highly persistent latent variable series. The novel approach of ancillarity-sufficiency interweaving (Yu and Meng, 2011) has recently been shown to aid in overcoming these issues for a broad class of multilevel models. This package provides software for ``combining best of different worlds'' which allows for inference for parameter constellations that have previously been infeasible to estimate without the need to select a particular parameterization beforehand.
}
\note{
This package is currently in active development. Your comments, suggestions and requests are warmly welcome!
}
\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}
\references{
Kastner, G. and Frühwirth-Schnatter, S. (2013). Ancillarity-sufficiency interweaving strategy (ASIS) for boosting MCMC estimation of stochastic volatility models. \emph{Computational Statistics and Data Analysis}, \url{http://dx.doi.org/10.1016/j.csda.2013.01.002}.

Yu, Y. and Meng, X.-L. (2011). To Center or Not to Center: That is Not the Question---An Ancillarity-Suffiency Interweaving Strategy (ASIS) for Boosting MCMC Efficiency. \emph{Journal of Computational and Graphical Statistics}, \bold{20}, 531--570.
}
\keyword{package}
\keyword{models}
\keyword{ts}

\examples{
## Simulate a highly persistent SV process 
sim <- svsim(500, mu = -10, phi = 0.99, sigma = 0.2)

## Obtain 5000 draws from the sampler (that's too little!)
draws <- svsample(sim$y, draws = 5000, burnin = 100, priormu = c(-10, 1),
                  priorphi = c(20, 1.2), priorsigma = 0.2)

## Predict 50 days ahead
fore <- predict(draws, 50)

## plot the results
plot(draws, forecast = fore)
}
