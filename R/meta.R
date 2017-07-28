#' Simple Meta-Analysis: Differential Expression
#' 
#' Perform differential expression analysis (using 
#' \code{\link[stats]{kruskal.test}}) over multiple replicate experiments, 
#' followed by Stouffer p-value combination and FDR analysis.
#' 
#' @param x an n by m numeric matrix, where m = num of samples, n = num of 
#'   molecules. Numerical values representing the expression level of the 
#'   molecule.
#' @param r factor. Sample replicate identity.
#' @param g factor. Sample cluster identity, passed to 
#'   \code{\link[stats]{kruskal.test}}.
#' @param var_thresh numeric. Molecules are only analyzed if their variance is 
#'   greater than this threshold in all replicate groups.
#'   
#' @export
#' 
#' @return a list with the following elements: \itemize{ \item{meta_pvals}{ 
#'   Stouffer p-values.} \item{FDR}{ FDR q-values.} \item{kruskal_pvals}{ 
#'   matrix of Kruskal test p-values used in meta-analysis.} 
#'   \item{is_replicated}{ logical. TRUE if molecule passes variance filter.} }
#'   
#' @importFrom stats kruskal.test var p.adjust
#' @importFrom metap sumz
#'   
#' @examples
#' library("SummarizedExperiment")
#' data("fluidigm",package = "scRNAseq")
#'
#' x = assay(fluidigm)[sample(x = seq_len(nrow(fluidigm)),size = 1000),]
#' r = factor(fluidigm$Coverage_Type)
#' g = factor(fluidigm$Biological_Condition)
#' kruskal_meta_obj = kruskalMeta(x,r,g)
#'
#' # Compare significance of 2 coverage types and color by reproducibility
#' plot(
#'   -log10(kruskal_meta_obj$kruskal_pvals)[kruskal_meta_obj$is_replicated, ],
#'   main = "-log10(p-value) by Coverage Type",
#'   xlim = c(0, max(
#'     -log10(kruskal_meta_obj$kruskal_pvals), na.rm = TRUE
#'   )),
#'   ylim = c(0, max(
#'     -log10(kruskal_meta_obj$kruskal_pvals), na.rm = TRUE
#'   )),
#'   col = 1 +
#'     (kruskal_meta_obj$FDR[kruskal_meta_obj$is_replicated] < 0.01),
#'   pch = 16
#' )
#' abline(0, 1, lty = 2)
#'

kruskalMeta = function(x,
                       r,
                       g,
                       var_thresh = 10 ^ -5) {
  is_var = rowSums(sapply(levels(r), function(p) {
    apply(x[, r == p], 1, var)
  }) > var_thresh) == nlevels(r)
  vx = x[is_var, ]
  
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
  
  meta_pvals <- FDR <- rep(NA, nrow(x))
  names(meta_pvals) <- names(FDR) <- rownames(x)
  meta_pvals[is_var] = apply(p_val_matrix[is_var, ],
                             1, function(x) {
                               metap::sumz(x)$p
                             })
  FDR[is_var] = stats::p.adjust(meta_pvals[is_var], method = "fdr")
  
  list(
    meta_pvals = meta_pvals,
    FDR = FDR,
    kruskal_pvals = p_val_matrix,
    is_replicated = is_replicated
  )
}

#' Simple Meta-Analysis: Signature Comparison
#' 
#' Compare distributions of signature values to a reference distribution (using
#' \code{\link[stats]{ks.test}}) over multiple replicate experiments, 
#' followed by Stouffer p-value combination and FDR analysis.
#' 
#' @param x an n by m numeric matrix, where m = num of samples, n = num of 
#'   signatures. Numerical values representing signature values.
#' @param x0 reference signature distribution.
#' @param r factor. Sample replicate identity.
#' @param var_thresh numeric. Molecules are only analyzed if their variance is 
#'   greater than this threshold in all replicate groups.
#' @param ... Additional arguments passed to \code{\link[stats]{ks.test}}).
#'   
#' @export
#' 
#' @return a list with the following elements: \itemize{ \item{meta_pvals}{ 
#'   Stouffer p-values.} \item{FDR}{ FDR q-values.} \item{ks_pvals}{ 
#'   matrix of KS-test p-values used in meta-analysis.} 
#'   \item{is_replicated}{ logical. TRUE if molecule passes variance filter.} }
#'   
#' @importFrom stats ks.test var p.adjust
#' @importFrom metap sumz
#'   
#' @examples
#' x = matrix(rnorm(1000),nrow = 10)
#' x0 = rnorm(1000)
#' r = factor(rep(1:2,each = 50))
#' ks_meta_obj = ksMeta(x = x,x0 = x0,r)
#'
#' # Compare significance of 2 coverage types and color by reproducibility
#' plot(-log10(ks_meta_obj$ks_pvals)[ks_meta_obj$is_replicated,],
#'      main = "-log10(p-value) by Replicate",
#'      xlim = c(0,max(-log10(ks_meta_obj$ks_pvals),na.rm = TRUE)),
#'      ylim = c(0,max(-log10(ks_meta_obj$ks_pvals),na.rm = TRUE)),
#'      col = 1 + (ks_meta_obj$FDR[ks_meta_obj$is_replicated] < 0.01),
#'      pch = 16)
#' abline(0,1, lty = 2)
#' 

ksMeta = function(x,
                  x0,
                  r,
                  var_thresh = 10 ^ -5,
                  ...) {
  is_var = rowSums(sapply(levels(r), function(p) {
    apply(x[, r == p], 1, var)
  }) > var_thresh) == nlevels(r)
  vx = x[is_var, ]
  
  p_val_matrix = matrix(
    NA,
    ncol = nlevels(r),
    nrow = nrow(x),
    dimnames = list(rownames(x), levels(r))
  )
  for (p in levels(r)) {
    is_p = r == p
    ks_list = apply(x[, is_p], 1, ks.test, ... , y = x0[is_p])
    p_val_matrix[, p] = unlist(lapply(ks_list,
                                      function(x) {
                                        x$p.value
                                      }))
  }
  
  is_replicated = !apply(is.na(p_val_matrix), 1, any)
  stopifnot(all(is_replicated == is_var))
  
  meta_pvals <- FDR <- rep(NA, nrow(x))
  names(meta_pvals) <- names(FDR) <- rownames(x)
  meta_pvals[is_var] = apply(p_val_matrix[is_var, ], 1, function(x) {
    metap::sumz(x)$p
  })
  FDR[is_var] = stats::p.adjust(meta_pvals[is_var], method = "fdr")
  
  list(
    meta_pvals = meta_pvals,
    FDR = FDR,
    ks_pvals = p_val_matrix,
    is_replicated = is_replicated
  )
}