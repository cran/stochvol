\encoding{UTF-8}
\name{.svsample}
\alias{.svsample}
\title{
Minimal overhead version of \code{\link{svsample}}.
}
\description{
 \code{.svsample} is a minimal overhead version of \code{\link{svsample}} with slightly different default arguments and a simplified return value structure. It is intended to be used mainly for one-step updates where speed is an issue, e.g., as a plug-in into other MCMC samplers. Note that absolutely no input checking is performed, thus this function is to be used with proper care!
}
\usage{
.svsample(y, draws = 1, burnin = 0, priormu = c(-10, 3), priorphi = c(5, 1.5),
	  priorsigma = 1, thinpara = 1, thinlatent = 1, thintime = 1,
	  quiet = TRUE, startpara, startlatent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
   numeric vector containing the data (usually log-returns), which must not contain zeroes.
}
  \item{draws}{
   single number greater or equal to 1, indicating the number of draws after burn-in (see below). Will be automatically coerced to integer. The defaults value is 1.
}
  \item{burnin}{
   single number greater or equal to 0, indicating the number of draws discarded as burn-in. Will be automatically coerced to integer. The default value is 0.
}
  \item{priormu}{
   numeric vector of length 2, indicating mean and standard deviation for the Gaussian prior distribution of the parameter \code{mu}, the level of the log-volatility. The default value is \code{c(-10, 3)}, which constitutes a rather vague prior for many common exchange rate datasets, stock returns and the like.
}
  \item{priorphi}{
   numeric vector of length 2, indicating the shape parameters for the Beta prior distributoin of the transformed parameter \code{(phi+1)/2}, where \code{phi} denotes the persistence of the log-volatility. The default value is \code{c(5, 1.5)}, which constitutes a prior that puts some belief in a persistent log-volatility but also encompasses the region where \code{phi} is around 0.
}
  \item{priorsigma}{
   single positive real number, which stands for the scaling of the transformed parameter \code{sigma^2}, where \code{sigma} denotes the volatility of log-volatility. More precisely, \code{sigma^2 ~ priorsigma * chisq(df = 1)}. The default value is \code{1}, which constitutes a reasonably vague prior for many common exchange rate datasets, stock returns and the like.
}
  \item{thinpara}{
   :single number greater or equal to 1, coercible to integer. Every \code{thinpara}th parameter draw is kept and returned. The default value is 1, corresponding to no thinning of the parameter draws -- every draw is stored.
}
  \item{thinlatent}{
   single number greater or equal to 1, coercible to integer. Every \code{thinlatent}th latent variable draw is kept and returned. The default value is 1, corresponding to no thinning of the latent variable draws, i.e. every draw is kept.
  }
  \item{thintime}{
   single number greater or equal to 1, coercible to integer. If \code{thintime} is different from 1, only every \code{thintime}th latent log-volatility is being monitored. If, e.g., \code{thintime = 3}, the latent log-volatilities \code{h_1,h_4,h_7,...} will be kept. The default value is 1, meaning that all latent variables \code{h_1,h_2,h_3,...} are stored.
}
  \item{quiet}{
   logical value indicating whether the progress bar and other informative output during sampling should be omitted. The default value is \code{TRUE}, implying non-verbose output.
}
  \item{startpara}{
   \emph{compulsory} named list, containing the starting values for the parameter draws. \code{startpara} must contain three elements named \code{mu}, \code{phi}, and \code{sigma}, where \code{mu} is an arbitrary numerical value, \code{phi} is a real number between \code{-1} and \code{1}, and \code{sigma} is a positive real number.
}
  \item{startlatent}{
   \emph{compulsory} vector of length \code{length(x$y)}, containing the starting values for the latent log-volatility draws.
}
}

\value{
A list with three components:
\item{para}{\code{3} times \code{draws} matrix containing the parameter draws.}
\item{latent}{\code{length(y)} times \code{draws} matrix containing draws of the latent variables \code{h_1, \dots, h_T}}.
\item{latent0}{Vector of length \code{draws} containing the draw(s) of the initial latent variable \code{h_0}.}
}
\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}

\details{
As opposed to the ordinary \code{\link{svsample}}, the default values differ for \code{draws}, \code{burnin}, and \code{quiet}. Note that currently neither \code{expert} nor \code{\dots} arguments are provided.
}

\note{
An example will soon be provided in a vignette.
}

\section{Warning}{
Expert use only! For most application, the use of \code{\link{svsample}} is recommended.
}
\seealso{
\code{\link{svsample}}
}
\keyword{models}
\keyword{ts}
