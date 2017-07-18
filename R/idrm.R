#' Log density of multivariate Gaussian distribution with symmetric
#' marginals
#'
#' Compute the log-density for 3-parameter multivariate Gaussian distribution.
#'
#' @param z data matrix
#' @param mu scalar mean
#' @param sigma	scalar standard deviation (diagonal covariance)
#' @param rho	scalar correlation coefficient (off-diagonal correlation)
#'
#' @export
#' @keywords internal
#'
#' @return Log density of multivariate Gaussian distribution.
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @examples
#' z.1 <- rnorm(500, 3, 1)
#' rho <- 0.8
#
#' ## The component with higher values is correlated with corr. coefficient=0.8
#' z.2 <- rnorm(500, 3 + 0.8*(z.1-3), (1-rho^2))
#' mu <- 3
#' sigma <- 1
#' den.z <- d.mnormal(cbind(z.1,z.2), mu, sigma, rho)
#'
#' den.z

d.mnormal = function(z, mu, sigma, rho)
{
  loglik <- mvtnorm::dmvnorm(
    x = z,
    mean = rep(mu, times = ncol(z)),
    sigma = (sigma ^ 2) * ((1 - rho) * diag(ncol(z)) + rho *
                             matrix(
                               1, nrow = ncol(z), ncol = ncol(z)
                             )),
    log = TRUE
  )
  return(loglik)
}

#' Compute log-likelihood of parameterized multivariate 2-component Gaussian
#' mixture models
#'
#' Compute the log-likelihood for parameterized multivariate 2-component
#' Gaussian mixture models with (1-p)N(0, 1, 0) + pN(mu, sigma, rho).
#'
#' @param z data matrix
#' @param mu scalar mean for the reproducible component.
#' @param sigma	scalar standard deviation (diagonal covariance) of the
#'   reproducible component.
#' @param rho	scalar correlation coefficient (off-diagonal correlation) of the
#'   reproducible component.
#' @param p	mixing proportion of the reproducible component.
#'
#' @export
#' @keywords internal
#'
#' @return Log-likelihood of the multivariate 2-component Gaussian mixture
#'   models (1-p)N(0, 1, 0) + pN(mu, sigma, rho).
#'
#' @examples
#' z.1 <- c(rnorm(500, 0, 1), rnorm(500, 3, 1))
#' rho <- 0.8
#'
#' ## The component with higher values is correlated with corr. coefficient=0.8
#' z.2 <- c(rnorm(500, 0, 1), rnorm(500, 3 + 0.8*(z.1[501:1000]-3), (1-rho^2)))
#'
#' ## Starting values
#' mu <- 3
#' sigma <- 1
#' rho <- 0.85
#' p <- 0.55
#'
#' ## The function is currently defined as
#' loglik <- loglik.2mnormal(cbind(z.1,z.2), mu, sigma, rho, p)
#'
#' loglik

loglik.2mnormal = function(z, mu, sigma, rho, p)
{
  l.m <-
    sum(d.mnormal(z, 0, 1, 0) + log(p * exp(
      d.mnormal(z, mu, sigma, rho) - d.mnormal(z, 0, 1, 0)
    ) +
      (1 - p)))
  return(l.m)
}

#' E-step for parameterized multivariate 2-component Gaussian mixture
#' models
#'
#' Expectation step in the EM algorithm for parameterized multivariate
#' 2-component Gaussian mixture models with (1-p)N(0, 1, 0) + pN(mu, sigma,
#' rho).
#'
#' @param z data matrix
#' @param mu scalar mean for the reproducible component.
#' @param sigma	scalar standard deviation (diagonal covariance) of the
#'   reproducible component.
#' @param rho	scalar correlation coefficient (off-diagonal correlation) of the
#'   reproducible component.
#' @param p	mixing proportion of the reproducible component.
#'
#' @export
#' @keywords internal
#'
#' @return e.z a numeric vector, where each entry represents the estimated
#'   expected conditional probability that an observation is in the
#'   reproducible component.
#'
#' @examples
#' z.1 <- c(rnorm(500, 0, 1), rnorm(500, 3, 1))
#' rho <- 0.8
#'
#' ## The component with higher values is correlated with corr. coefficient=0.8
#' z.2 <- c(rnorm(500, 0, 1), rnorm(500, 3 + 0.8*(z.1[501:1000]-3), (1-rho^2)))
#'
#' ## Starting values
#' mu0 <- 3
#' sigma0 <- 1
#' rho0 <- 0.85
#' p0 <- 0.55
#'
#' e.z <- e.step.2mnormal(cbind(z.1,z.2), mu0, sigma0, rho0, p0)

e.step.2mnormal = function (z, mu, sigma, rho, p)
{
  e.z <-
    p / ((1 - p) * exp(d.mnormal(z, 0, 1, 0) -
                         d.mnormal(z, mu, sigma, rho)) + p)
  invisible(e.z)
}

#' M-step for parameterized multivariate 2-component Gaussian mixture
#' models
#'
#' Maximization step in the EM algorithm for parameterized multivariate
#' 2-component Gaussian mixture models with (1-p)N(0, 1, 0) + pN(mu, sigma,
#' rho).
#'
#' @param z data matrix
#' @param e.z a vector of expected conditional probability that the
#'   corresponding observation (Z row) is reproducible.
#'
#' @export
#' @keywords internal
#'
#' @return a list with the following elements: \itemize{ \item{p}{the estimated
#'   mixing proportion of the reproducible component.} \item{mu}{the estimated
#'   scalar mean for the reproducible component.} \item{sigma}{the estimated
#'   scalar standard deviation (diagonal covariance) of the reproducible
#'   component.} \item{rho}{the estimated scalar correlation coefficient
#'   (off-diagonal correlation) of the reproducible component.} }
#'
#' @examples
#' z.1 <- c(rnorm(500, 0, 1), rnorm(500, 3, 1))
#' rho <- 0.8
#'
#' ##The component with higher values is correlated with corr. coefficient=0.8
#' z.2 <- c(rnorm(500, 0, 1),
#'          rnorm(500, 3 + 0.8 * (z.1[501:1000] - 3), (1 - rho ^ 2)))
#'
#' e.z <- c(rep(0, 500) + abs(rnorm(500, 0, 0.05)),
#'          rep(1, 500) - abs(rnorm(500, 0, 0.05)))
#' para <- m.step.2mnormal(cbind(z.1, z.2), e.z)
#' para

m.step.2mnormal = function(z, e.z)
{
  p <- mean(e.z)
  mu <- sum(rowMeans(z) * e.z) / sum(e.z)
  centered_z = (z - mu) * sqrt(e.z)
  C = t(centered_z) %*% centered_z
  sigma <- sqrt(sum(diag(C)) / (sum(e.z) * ncol(z)))
  rho <- (sum(C) - sum(diag(C))) / (sum(diag(C)) * (ncol(z) - 1))
  return(list(
    p = p,
    mu = mu,
    sigma = sigma,
    rho = rho
  ))
}

#' Estimate the irreproducible discovery rate using the copula mixture
#' model
#'
#' Fit a multivariate Gaussian copula mixture model.
#'
#' @param x an n by m numeric matrix, where m = num of replicates, n = num of
#'   observations. Numerical values representing the significance of the
#'   observations. Note that significant signals are expected to have large
#'   values of x. In case that smaller values represent higher significance
#'   (e.g. p-value), a monotonic transformation needs to be applied to reverse
#'   the order before using this function, for example, -log(p-value).
#' @param mu a starting value for the scalar mean for the reproducible
#'   component.
#' @param sigma	a starting value for the scalar standard deviation (diagonal
#'   covariance) of the reproducible component.
#' @param rho	a starting value for the scalar correlation coefficient
#'   (off-diagonal correlation) of the reproducible component.
#' @param p	a starting value for the proportion of the reproducible component.
#' @param eps	Stopping criterion. Iterations stop when the increment of
#'   log-likelihood is < eps*log-likelihood, Default=0.001.
#' @param max.ite	Maximum number of iterations. Default=30.
#'
#' @export
#'
#' @importFrom stats ecdf
#'
#' @return a list with the following elements: \itemize{ \item{para}{ estimated
#'   parameters: p, rho, mu, sigma.} \item{idr}{ a numeric vector of the local
#'   idr for each observation (i.e. estimated conditional probablility for each
#'   observation to belong to the irreproducible component.} \item{IDR}{ a
#'   numerical vector of the expected irreproducible discovery rate for
#'   observations that are as irreproducible or more irreproducible than the
#'   given observations.} \item{loglik}{ log-likelihood at the end of
#'   iterations.} \item{loglik.trace}{ trajectory of log-likelihood.} }
#'
#' @examples
#' data("simu.idr",package = "idr")
#' # simu.idr$x and simu.idr$y are p-values
#' # Transfer them such that large values represent significant ones
#' x <- cbind(-simu.idr$x, -simu.idr$y)
#'
#' mu <- 2.6
#' sigma <- 1.3
#' rho <- 0.8
#' p <- 0.7
#'
#' idr.out <- est.IDRm(x, mu, sigma, rho, p, eps=0.001, max.ite=20)
#'
#' names(idr.out)

est.IDRm = function(x,
                    mu,
                    sigma,
                    rho,
                    p,
                    eps = 0.001,
                    max.ite = 30)
{
  conv <- function(old, new){
    if(is.infinite(new)){
      stop("likelihood is infinite: convergence condition invalid")
    }
    return(abs(new - old) < eps * (1 + abs(new)))
  }
  
  # Steps 1-2: Compute cdfs and rescale
  x.cdf.funcs = apply(x, 2, ecdf)
  afactor <- nrow(x) / (nrow(x) + 1)
  x.cdf = sapply(seq_len(ncol(x)), function(s) {
    x.cdf.funcs[[s]](x[, s]) * afactor
  })
  
  # Step 3: Initialize parameters
  para <- list()
  para$mu <- mu
  para$sigma <- sigma
  para$rho <- rho
  para$p <- p
  
  j <- 1
  to.run <- TRUE
  loglik.trace <- c()
  loglik.inner.trace <- c()
  z <-
    apply(
      x.cdf,
      2,
      idr::get.pseudo.mix,
      mu = para$mu,
      sigma = para$sigma,
      rho = para$rho,
      p = para$p
    )
  while (to.run) {
    i <- 1
    while (to.run) {
      e.z <- e.step.2mnormal(z, para$mu, para$sigma,
                             para$rho, para$p)
      if(all(e.z == 0)){
        stop("all memberships numerically zero. try new initial params")
      }
      para <- m.step.2mnormal(z, e.z)
      if (i > 1)
        l.old <- l.new
      l.new <- loglik.2mnormal(z, para$mu, para$sigma,
                               para$rho, para$p)
      loglik.inner.trace[i] <- l.new
      if (i > 1) {
          to.run <- !conv(loglik.inner.trace[i - 1], loglik.inner.trace[i])
      }
      i <- i + 1
    }
    z <-
      apply(
        x.cdf,
        2,
        idr::get.pseudo.mix,
        mu = para$mu,
        sigma = para$sigma,
        rho = para$rho,
        p = para$p
      )
    if (j > 1)
      l.old.outer <- l.new.outer
    l.new.outer <- loglik.2mnormal(z, para$mu, para$sigma,
                                   para$rho, para$p)
    loglik.trace[j] <- l.new.outer
    if (j == 1)
      to.run <- TRUE
    else {
      if (j > max.ite)
        to.run <- FALSE
      else
        to.run <- !conv(l.old.outer, l.new.outer)
    }
    j <- j + 1
  }
  idr <- 1 - e.z
  o <- order(idr)
  idr.o <- idr[o]
  idr.rank <- rank(idr.o, ties.method = "max")
  top.mean <- function(index, x) {
    mean(x[1:index])
  }
  IDR.o <- sapply(idr.rank, top.mean, idr.o)
  IDR <- idr
  IDR[o] <- IDR.o
  return(list(
    para = list(
      p = para$p,
      rho = para$rho,
      mu = para$mu,
      sigma = para$sigma
    ),
    loglik = l.new,
    loglik.trace = loglik.trace,
    idr = 1 - e.z,
    IDR = IDR
  ))
}

#' Irreproducible Discovery Rate analysis: Differential Expression
#'
#' Perform differential expression analysis (using
#' \code{\link[stats]{kruskal.test}}) over multiple replicate experiments,
#' followed by irreproducible discovery rate (IDR) analysis.
#'
#' @param x an n by m numeric matrix, where m = num of samples, n = num of
#'   molecules. Numerical values representing the expression level of the
#'   molecule.
#' @param r factor. Sample replicate identity.
#' @param g factor. Sample cluster identity, passed to
#'   \code{\link[stats]{kruskal.test}}.
#' @param var_thresh numeric. Molecules are only analyzed if their variance is
#'   greater than this threshold in all replicate groups.
#' @param idr_mu a starting value for idr mean (see
#'   \code{\link[scrap]{est.IDRm}}).
#' @param idr_sigma	a starting value for the idr standard deviation (see
#'   \code{\link[scrap]{est.IDRm}}).
#' @param idr_rho	a starting value for the idr correlation coefficient (see
#'   \code{\link[scrap]{est.IDRm}}).
#' @param idr_p	a starting value for the proportion of the reproducible tests
#'   (see \code{\link[scrap]{est.IDRm}}).
#'
#' @export
#'
#' @return a list with the following elements: \itemize{ \item{idr}{ list.
#'   \code{\link[scrap]{est.IDRm}} output.}\item{kruskal_pvals}{ matrix of
#'   Kruskal test p-values used in IDR analysis.} \item{is_replicated}{
#'   logical. TRUE if molecule passes variance filter.} }
#'
#' @importFrom stats kruskal.test var
#'
#' @examples
#' library("SummarizedExperiment")
#' data("fluidigm",package = "scRNAseq")
#'
#' x = assay(fluidigm)[sample(x = seq_len(nrow(fluidigm)),size = 1000),]
#' r = factor(fluidigm$Coverage_Type)
#' g = factor(fluidigm$Biological_Condition)
#' kruskal_idrm_obj = kruskalIDRm(x,r,g)
#'
#' # Compare significance of 2 coverage types and color by reproducibility
#' plot(-log10(kruskal_idrm_obj$kruskal_pvals)[kruskal_idrm_obj$is_replicated,],
#'      main = "-log10(p-value) by Coverage Type",
#'      xlim = c(0,max(-log10(kruskal_idrm_obj$kruskal_pvals),na.rm = TRUE)),
#'      ylim = c(0,max(-log10(kruskal_idrm_obj$kruskal_pvals),na.rm = TRUE)),
#'      col = 1 + (kruskal_idrm_obj$idr$IDR < 0.01), pch = 16)
#' abline(0,1, lty = 2)
#'

kruskalIDRm = function(x,
                       r,
                       g,
                       var_thresh = 10 ^ -5,
                       idr_mu = 2,
                       idr_sigma = 2,
                       idr_rho = .5,
                       idr_p = 0.01) {
  is_var = rowSums(sapply(levels(r), function(p) {
    apply(x[, r == p], 1, var)
  }) > var_thresh) == nlevels(r)
  vx = x[is_var,]
  
  p_val_matrix = matrix(
    NA,
    ncol = nlevels(r),
    nrow = nrow(x),
    dimnames = list(rownames(x), levels(r))
  )
  for (p in levels(r)) {
    is_p = r == p
    kruskal_list = apply(vx[, is_p], 1,
                         kruskal.test, g = g[is_p])
    p_val_matrix[, p][is_var]  = unlist(lapply(kruskal_list, function(x) {
      x$p.value
    }))
  }
  
  is_replicated = !apply(is.na(p_val_matrix), 1, any)
  stopifnot(all(is_replicated == is_var))
  idr_objm = est.IDRm(
    -log10(p_val_matrix[is_replicated, ]),
    mu = idr_mu,
    sigma = idr_sigma,
    rho = idr_rho,
    p = idr_p
  )
  list(idr = idr_objm,
       kruskal_pvals = p_val_matrix,
       is_replicated = is_replicated)
}