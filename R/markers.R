#' Reproducible Marker Prediction
#' 
#' Wrapper function for differential expression (de) analysis and reproducible
#' module (module) analysis.
#' 
#' @param x an n by m numeric matrix, where m = num of samples, n = num of 
#'   molecules. Numerical values representing the expression level of the 
#'   molecule.
#' @param r factor. Sample replicate identity.
#' @param g factor. Sample cluster identity.
#' @param external_markers character. If not NULL, only return markers if they
#'   are included in this argument.
#' @param de_idr_thresh numeric. IDR threshold for calling differential
#'   expression (default 0.01).
#' @param module_bonf_thresh numeric. Bonferroni-adjusted threshold for calling
#'   reproducible module genes (default 0.01).
#' @param var_thresh numeric. Molecules are only analyzed if their variance is 
#'   greater than this threshold in all replicate groups.
#' @param mad_constant scale factor passed to \code{\link[stats]{mad}} in 
#'   scaling z-transformed correlation coefficients.
#' @param pair_pthresh numeric. Molecule pairs are called "reproducibly 
#'   correlated" if, across all replicate groups, their z-value corresponds to 
#'   a two-sided p-value below this threshold.
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
#' @return character. Genes in the intersection of de, modules, and (optional)
#'   external marker list.
#'   
#' @importFrom igraph graph.adjacency
#' @importFrom stats p.adjust
#'   
#' @examples
#' library("SummarizedExperiment")
#' data("fluidigm",package = "scRNAseq")
#'
#' x = assay(fluidigm)[sample(x = seq_len(nrow(fluidigm)),size = 1000),]
#' r = factor(fluidigm$Coverage_Type)
#' g = factor(fluidigm$Biological_Condition)
#' getMarkers(x,r,g,de_idr_thresh = 0.01)
#'

getMarkers = function(x,
                      r,
                      g,
                      external_markers = NULL,
                      de_idr_thresh = 0.01,
                      module_bonf_thresh = 0.01,
                      var_thresh = 10 ^ -5,
                      mad_constant = 1.4826,
                      pair_pthresh = 0.1,
                      idr_mu = 2,
                      idr_sigma = 2,
                      idr_rho = .5,
                      idr_p = 0.01) {
  
  # DE
  kruskal_idrm_obj = kruskalIDRm(
    x = x,
    r = r,
    g = g,
    var_thresh = var_thresh,
    idr_mu = idr_mu,
    idr_sigma = idr_sigma,
    idr_rho = idr_rho,
    idr_p = idr_p
  )
  
  all_names = rownames(x)[kruskal_idrm_obj$is_replicated]
  de_names = all_names[kruskal_idrm_obj$idr$IDR < de_idr_thresh]
  
  # Modules
  A = get.repro.thresh.adjacency(
    x = x,
    r = r,
    mad_constant = mad_constant,
    var_thresh = var_thresh,
    pair_pthresh = pair_pthresh
  )
  pdegree = pzipdegree(graph.adjacency(A, mode = "undirected"))
  module_names = names(which(p.adjust(pdegree,
                                      method = "bonf") < module_bonf_thresh))
  
  # Markers
  marker_names = sort(intersect(module_names, de_names))
  
  if (!is.null(external_markers)) {
    marker_names = sort(intersect(external_markers, marker_names))
  }
  
  marker_names
}