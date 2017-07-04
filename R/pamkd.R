#' Partitioning around the medoids with PCA dimension reduction
#' 
#' This calls the function \code{\link[fpc]{pamk}} to perform a partitioning 
#' around medoids clustering with the number of clusters estimated by optimum 
#' average silhouette width. Silhouette widths are computed on a Euclidean 
#' metric over the first \code{d} principal components, where \code{d} is 
#' varied over a range for each run. The optimal run is defined as the run with
#' maximal number of clusters over the range, using average silhouette width to
#' break ties.
#' 
#' @param x a numeric or complex matrix (or data frame) which provides the data
#'   for the principal components analysis.
#' @param maxk numeric. The maximum number of clusters to be passed to 
#'   \code{\link[fpc]{pamk}}
#' @param maxd numeric. The maximum index of prinicipal components, \code{d}, 
#'   to be inluded in the distance metric.
#' @param to_log logical. If TRUE (default FALSE), apply 
#'   \code{\link[base]{log1p}} transformation prior to PCA (via 
#'   \code{\link[rARPACK]{svds}}).
#' @param pca_center logical. If TRUE (default), apply row-centering prior to 
#'   \code{\link[rARPACK]{svds}}.
#' @param pca_scale logical. If TRUE (default), apply row-scaling prior to 
#'   \code{\link[rARPACK]{svds}}.
#' @param verbose logical. If TRUE (default FALSE), print informative messages.
#'   
#' @export
#' 
#' @return a list with the following elements: \itemize{ \item{pamobject}{ The 
#'   output of the optimal run of the \code{\link[cluster]{pam}}-function.} 
#'   \item{nc}{ the optimal number of clusters.} \item{nd}{ the optimal number
#'   of principal components.} \item{Y}{ dimension reduced data with optimal 
#'   dimension.} }
#'   
#' @importFrom fpc pamk
#' @importFrom rARPACK svds
#' @importFrom cluster pam
#' @importFrom methods as
#' @importFrom stats dist
#'   
#' @examples
#' x <- cbind(matrix(rpois(10000, lambda = 3), ncol = 100),
#'            matrix(rpois(10000, lambda = 4), ncol = 100))
#' pamkd_out = pamkd(x, to_log = TRUE, pca_scale = FALSE)
#' 
#' # Show the objects in the 2D tsne representation
#' # tsne_obj = Rtsne(pamkd_out$Y, pca = FALSE) # Run TSNE
#' # plot(tsne_obj$Y,
#' #     col = pamkd_out$pamobject$clustering,
#' #     pch = 16)

pamkd = function(x,
                 maxk = min(5, min(dim(x)) - 1),
                 maxd = min(5, min(dim(x)) - 1),
                 to_log = FALSE,
                 pca_center = TRUE,
                 pca_scale = TRUE,
                 verbose = FALSE) {
  if (verbose) {
    print(paste0("maxk = ", maxk))
    print(paste0("maxd = ", maxd))
  }
  if (to_log) {
    x = log1p(x)
  }
  x = apply(x, 1, scale, scale = pca_scale, center = pca_center)
  svds_obj = rARPACK::svds(
    methods::as(x, "dgCMatrix"),
    k = maxd,
    nu = maxd,
    nv = 0
  )
  pc_mat = t(t(svds_obj$u) * svds_obj$d)
  nc_all = rep(-1, maxd)
  crit_all = rep(-1, maxd)
  totdist = stats::dist(pc_mat[, 1]) ^ 2
  for (d in 2:maxd) {
    totdist = totdist + stats::dist(pc_mat[, d]) ^ 2
    pamkobj = fpc::pamk(sqrt(totdist), krange = 1:maxk)
    nc_all[d] = pamkobj$nc
    crit_all[d] = pamkobj$crit[pamkobj$nc]
    if (verbose) {
      print(paste0("Optimized nc for ", d, " dimensions"))
    }
  }
  
  nc = max(nc_all)
  crit_nc = crit_all[which(nc_all ==  nc)]
  nd = which(nc_all ==  nc)[crit_nc == max(crit_nc)]
  
  list(
    pamobject = cluster::pam(stats::dist(pc_mat[, 1:nd]), k = nc),
    nc = nc,
    nd = nd,
    Y = pc_mat[, 1:nd]
  )
  
}