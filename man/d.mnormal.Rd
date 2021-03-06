% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idrm.R
\name{d.mnormal}
\alias{d.mnormal}
\title{Log density of multivariate Gaussian distribution with symmetric
marginals}
\usage{
d.mnormal(z, mu, sigma, rho)
}
\arguments{
\item{z}{data matrix}

\item{mu}{scalar mean}

\item{sigma}{scalar standard deviation (diagonal covariance)}

\item{rho}{scalar correlation coefficient (off-diagonal correlation)}
}
\value{
Log density of multivariate Gaussian distribution.
}
\description{
Compute the log-density for 3-parameter multivariate Gaussian distribution.
Code adapted from \code{\link[idr]{d.binormal}}.
}
\examples{
z.1 <- rnorm(500, 3, 1)
rho <- 0.8
## The component with higher values is correlated with corr. coefficient=0.8
z.2 <- rnorm(500, 3 + 0.8*(z.1-3), (1-rho^2))
mu <- 3
sigma <- 1
den.z <- d.mnormal(cbind(z.1,z.2), mu, sigma, rho)

den.z
}
\keyword{internal}
