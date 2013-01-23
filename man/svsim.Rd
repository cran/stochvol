\name{svsim}
\alias{svsim}
\title{Simulating a Stochastic Volatility Process}

\description{
 \code{svsim} is used to produce realizations of a stochastic volatility (SV) process.
}

\usage{
svsim(len, mu = -10, phi = 0.98, sigma = 0.2)
}

\arguments{
  \item{len}{
   length of the simulated time series.
  }
  \item{mu}{
   level of the latent log-volatility AR(1) process. The defaults value is -10.
  }
  \item{phi}{
   persistence of the latent log-volatility AR(1) process. The default value is 0.98.
  }
  \item{sigma}{
   volatility of the latent log-volatility AR(1) process. The defaults value is 0.2.
  }
}

\details{
 This function draws an initial log-volatily \code{h_0} from the stationary distribution of the AR(1) process and iteratively generates \code{h_1,...,h_T}. Finally, the ``log-returns'' \code{y ~ N(0, exp(h/2))} are simulated.
}

\value{
The output is a list object of class \code{svsim} containing
\item{y}{a vector of length \code{len} containing the simulated data, usually interpreted as ``log-returns''.}
\item{vol}{a vector of length \code{len} containing the simulated instantaneous volatilities \code{exp(h_t/2)}.}
\item{vol0}{the initial volatility \code{exp(h_0/2)}, drawn from the stationary distribution of the latent AR(1) process.}
\item{para}{a named list with three elements \code{mu}, \code{phi}, and \code{sigma}, containing the corresponding arguments.}

To display the output use \code{print}, \code{summary} and \code{plot}. The \code{print} method simply prints the content of the object in a moderately formatted manner. The \code{summary} method provides some summary statistics (in \%), and the \code{plot} method plots the the simulated 'log-returns' \code{y} along with the corresponding volatilites \code{vol}.
}
\author{Gregor Kastner \email{gregor.kastner@wu.ac.at}
}

\seealso{\code{\link{svsample}}
}

\examples{
## Simulate a highly persistent SV process of length 500
sim <- svsim(500, phi = 0.99, sigma = 0.1)

print(sim)
summary(sim)
plot(sim)
}
\keyword{datagen}
\keyword{ts}
