% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idrm.R
\name{m.step.2mnormal}
\alias{m.step.2mnormal}
\title{M-step for parameterized multivariate 2-component Gaussian mixture
models}
\usage{
m.step.2mnormal(z, e.z)
}
\arguments{
\item{z}{data matrix}

\item{e.z}{a vector of expected conditional probability that the
corresponding observation (Z row) is reproducible.}
}
\value{
a list with the following elements: \itemize{ \item{p}{the estimated
  mixing proportion of the reproducible component.} \item{mu}{the estimated
  scalar mean for the reproducible component.} \item{sigma}{the estimated
  scalar standard deviation (diagonal covariance) of the reproducible
  component.} \item{rho}{the estimated scalar correlation coefficient
  (off-diagonal correlation) of the reproducible component.} }
}
\description{
Maximization step in the EM algorithm for parameterized multivariate
2-component Gaussian mixture models with (1-p)N(0, 1, 0) + pN(mu, sigma,
rho).
Code adapted from \code{\link[idr]{m.step.2normal}}.
}
\examples{
z.1 <- c(rnorm(500, 0, 1), rnorm(500, 3, 1))
rho <- 0.8

##The component with higher values is correlated with corr. coefficient=0.8
z.2 <- c(rnorm(500, 0, 1),
         rnorm(500, 3 + 0.8 * (z.1[501:1000] - 3), (1 - rho ^ 2)))

e.z <- c(rep(0, 500) + abs(rnorm(500, 0, 0.05)),
         rep(1, 500) - abs(rnorm(500, 0, 0.05)))
para <- m.step.2mnormal(cbind(z.1, z.2), e.z)
para
}
\keyword{internal}
