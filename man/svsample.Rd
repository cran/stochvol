\encoding{UTF-8}
\name{svsample}
\alias{svsample}
\title{
Markov Chain Monte Carlo (MCMC) Sampling for the Stochastic Volatility (SV) Model
}
\description{
\code{svsample} simulates from the joint posterior distribution of the SV parameters \code{mu}, \code{phi}, and \code{sigma}, along with the latent log-volatilities \code{h_0,...,h_n} and returns the MCMC draws.
}
\usage{
svsample(y, draws = 10000, burnin = 1000, priormu = c(0, 100),
         priorphi = c(5, 1.5), priorsigma = 1, thinpara = 1,
         thinlatent = 1, thintime = 1, quiet = FALSE, startpara,
         startlatent, expert, ...)
}
\arguments{
  \item{y}{
   numeric vector containing the data (usually log-returns), which must not contain zeroes. Alternatively, \code{y} can be an \code{svsim} object. In this case, the returns will be extracted and a warning is thrown.
}
  \item{draws}{
   single number greater or equal to 1, indicating the number of draws after burn-in (see below). Will be automatically coerced to integer. The defaults value is 10000.
}
  \item{burnin}{
   single number greater or equal to 0, indicating the number of draws discarded as burn-in. Will be automatically coerced to integer. The default value is 1000.
}
  \item{priormu}{
   numeric vector of length 2, indicating mean and standard deviation for the Gaussian prior distribution of the parameter \code{mu}, the level of the log-volatility. The default value is \code{c(0, 100)}, which constitutes a practically uninformative prior for common exchange rate datasets, stock returns and the like.
}
  \item{priorphi}{
   numeric vector of length 2, indicating the shape parameters for the Beta prior distribution of the transformed parameter \code{(phi+1)/2}, where \code{phi} denotes the persistence of the log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a prior that puts some belief in a persistent log-volatility but also encompasses the region where \code{phi} is around 0.
}
  \item{priorsigma}{
   single positive real number, which stands for the scaling of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma * chisq(df = 1)}. The default value is \code{1}, which constitutes a reasonably vague prior for many common exchange rate datasets, stock returns and the like.
}
  \item{thinpara}{
   single number greater or equal to 1, coercible to integer. Every \code{thinpara}th parameter draw is kept and returned. The default value is 1, corresponding to no thinning of the parameter draws i.e. every draw is stored.
}
  \item{thinlatent}{
   single number greater or equal to 1, coercible to integer. Every \code{thinlatent}th latent variable draw is kept and returned. The default value is 1, corresponding to no thinning of the latent variable draws, i.e. every draw is kept.
  }
  \item{thintime}{
   single number greater or equal to 1, coercible to integer. If \code{thintime} is different from 1, only every \code{thintime}th latent log-volatility is being monitored. If, e.g., \code{thintime = 3}, the latent log-volatilities \code{h_1,h_4,h_7,...} will be kept. The default value is 1, meaning that all latent variables \code{h_1,h_2,h_3,...} are stored.
}
  \item{quiet}{
   logical value indicating whether the progress bar and other informative output during sampling should be omitted. The default value is \code{FALSE}, implying verbose output.
}
  \item{startpara}{
   \emph{optional} named list, containing the starting values for the parameter draws. If supplied, \code{startpara} must contain three elements named \code{mu}, \code{phi}, and \code{sigma}, where \code{mu} is an arbitrary numerical value, \code{phi} is a real number between \code{-1} and \code{1}, and \code{sigma} is a positive real number. The default value is \code{list(mu = -10, phi = 0.9, sigma = 0.3)}.
}
  \item{startlatent}{
   \emph{optional} vector of length \code{length(x$y)}, containing the starting values for the latent log-volatility draws. The default value is \code{rep(-10, length(x$y))}.
}
  \item{expert}{
   \emph{optional} named list of expert parameters. For most applications, the default values probably work best. Interested users are referred to the literature provided in the References section. If \code{expert} is provided, it may contain the following named elements:
   
   \code{parameterization}: Character string equal to \code{"centered"}, \code{"noncentered"}, \code{"GIS_C"}, or \code{"GIS_NC"}. Defaults to \code{"GIS_C"}.

   \code{mhcontrol}: Single numeric value controlling the proposal density of a Metropolis-Hastings (MH) update step when sampling \code{sigma}. If \code{mhcontrol} is smaller than 0, an independence proposal will be used, while values greater than zero control the stepsize of a log-random-walk proposal. Defaults to \code{-1}.

   \code{gammaprior}: Single logical value indicating whether a Gamma prior for \code{sigma^2} should be used. If set to \code{FALSE}, an Inverse Gamma prior is employed. Defaults to \code{TRUE}.

   \code{truncnormal}: Single logical value indicating whether a truncated Gaussian distribution should be used as proposal for draws of \code{phi}. If set to \code{FALSE}, a regular Gaussian prior is employed and the draw is immediately discarded when values outside the unit ball happen to be drawn. Defaults to \code{FALSE}.

   \code{mhsteps}: Either \code{1}, \code{2}, or \code{3}. Indicates the number of blocks used for drawing from the posterior of the parameters. Defaults to \code{2}.

   \code{proposalvar4sigmaphi}: Single positive number indicating the conditional prior variance of \code{sigma*phi} in the ridge \emph{proposal} density for sampling \code{(mu, phi)}. Defaults to \code{10^8}.
   
   \code{proposalvar4sigmatheta}: Single positive number indicating the conditional prior variance of \code{sigma*theta} in the ridge \emph{proposal} density for sampling \code{(mu, phi)}. Defaults to \code{10^12}.
}
  \item{\dots}{
   Any extra arguments will be forwarded to \code{\link{updatesummary}}, controlling the type of statistics calculated for the posterior draws.
}
}
\details{
For details concerning the algorithm please see the paper by Kastner and Frühwirth-Schnatter (2013).
}
\value{
The value returned is a list object of class \code{svdraws} holding
\item{para}{\code{mcmc} object containing the \emph{parameter} draws from the posterior distribution.}
\item{latent}{\code{mcmc} object containing the \emph{latent instantaneous log-volatility} draws from the posterior distribution.}
\item{latent0}{\code{mcmc} object containing the \emph{latent initial log-volatility} draws from the posterior distribution.}
\item{y}{the argument \code{y}}
\item{runtime}{\code{proc_time} object containing the run time of the sampler.}
\item{priors}{\code{list} containing the parameter values of the prior distribution, i.e. the arguments \code{priormu}, \code{priorphi} and \code{priorsigma}.}
\item{thinning}{\code{list} containing the thinning parameters, i.e. the arguments \code{thinpara}, \code{thinlatent} and \code{thintime}.}
\item{summary}{\code{list} containing a collection of summary statistics of the posterior draws for \code{para}, \code{latent}, and \code{latent0}.}

To display the output, use \code{print}, \code{summary} and \code{plot}. The \code{print} method simply prints the posterior draws (which is very likely a lot of output); the \code{summary} method displays the summary statistics currently stored in the object; the \code{plot} method \code{\link{plot.svdraws}} gives a graphical overview of the posterior distribution by calling \code{\link{volplot}}, \code{\link{traceplot}} and \code{\link{densplot}} and displaying the results on a single page. 
}
\references{
Kastner, G. and Frühwirth-Schnatter, S. (2014). Ancillarity-sufficiency interweaving strategy (ASIS) for boosting MCMC estimation of stochastic volatility models. \emph{Computational Statistics & Data Analysis}, \bold{76}, 408--423, \url{http://dx.doi.org/10.1016/j.csda.2013.01.002}.
}
\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}
\note{
 If \code{y} contains zeros, you might want to consider de-meaning your returns.
}

\seealso{
\code{\link{svsim}}, \code{\link{updatesummary}}, \code{\link{predict.svdraws}}, \code{\link{plot.svdraws}}.
}
\examples{
## Simulate a short and highly persistent SV process 
sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)

## Obtain 5000 draws from the sampler (that's not a lot)
draws <- svsample(sim$y, draws = 5000, burnin = 100,
		  priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)

## Check out the results
summary(draws)
plot(draws)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{ts}
