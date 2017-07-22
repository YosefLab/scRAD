#' Reproducible Marker Prediction
#' 
#' Wrapper function for differential expression (de) analysis and reproducible
#' module (module) analysis.
#' 
#' @param x an n by m numeric matrix, where m = num of samples, n = num of 
#'   molecules. Numerical values representing the expression level of the 
#'   molecule.
#' @param r factor. Sample replicate identity.
#' @param g factor. Sample cluster identity. NA samples are ignored in de.
#' @param r_de factor. If specified, use this argument as de analysis sample
#'   replicate identity instead of \code{r}.
#' @param g_target character. If not null, only include genes up-regulated in
#'   this cluster.
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
#'   \code{\link[scider]{est.IDRm}}).
#' @param idr_sigma	a starting value for the idr standard deviation (see 
#'   \code{\link[scider]{est.IDRm}}).
#' @param idr_rho	a starting value for the idr correlation coefficient (see 
#'   \code{\link[scider]{est.IDRm}}).
#' @param idr_p	a starting value for the proportion of the reproducible tests 
#'   (see \code{\link[scider]{est.IDRm}}).
#' @param plot_venn	if TRUE, plot venn diagram of marker candidates.
#' @param ...	additional arguments passed to \code{\link[scider]{est.IDRm}}.
#'   
#' @export
#' 
#' @return character. Genes in the intersection of de, modules, and (optional)
#'   external marker list.
#'   
#' @importFrom igraph graph.adjacency
#' @importFrom stats p.adjust aggregate
#' @importFrom VennDiagram draw.triple.venn
#'   
#' @examples
#' library("SummarizedExperiment")
#' data("fluidigm",package = "scRNAseq")
#'
#' x = assay(fluidigm)[sample(x = seq_len(nrow(fluidigm)),size = 1000),]
#' r = factor(fluidigm$Coverage_Type)
#' g = factor(fluidigm$Biological_Condition)
#' getMarkers(x,r,g)
#'

getMarkers = function(x,
                      r,
                      g,
                      r_de = r,
                      g_target = NULL,
                      external_markers = NULL,
                      de_idr_thresh = 0.01,
                      module_bonf_thresh = 0.01,
                      var_thresh = 10 ^ -5,
                      mad_constant = 1.4826,
                      pair_pthresh = 0.01,
                      idr_mu = 2,
                      idr_sigma = 2,
                      idr_rho = .5,
                      idr_p = 0.01,
                      plot_venn = FALSE,
                      ...
                      ) {
  
  # DE
  is_compared = !is.na(g)
  kruskal_idrm_obj = kruskalIDRm(
    x = x[,is_compared],
    r = r_de[is_compared],
    g = g[is_compared],
    var_thresh = var_thresh,
    idr_mu = idr_mu,
    idr_sigma = idr_sigma,
    idr_rho = idr_rho,
    idr_p = idr_p,
    ...
  )
  
  all_names = rownames(x)[kruskal_idrm_obj$is_replicated]
  de_names = all_names[kruskal_idrm_obj$idr$IDR < de_idr_thresh]
  
  if(!is.null(g_target)){
    is_target = g == g_target
    agg_out = aggregate(t(x[,is_compared]),
                    by = list(is_target[is_compared]),
                    mean)
    up_names = rownames(x)[agg_out[,-1][2,] > agg_out[,-1][1,]]
    de_names = intersect(de_names,up_names)
  }
  
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
  
  
  # Venn
  if (plot_venn) {
    d1 = module_names
    d2 = external_markers
    d3 = de_names
    
    area1 = length(d1)
    area2 = length(d2)
    area3 = length(d3)
    n12 = length(intersect(d1, d2))
    n13 = length(intersect(d1, d3))
    n23 = length(intersect(d2, d3))
    n123 = length(intersect(d1, intersect(d2, d3)))
    
    VennDiagram::draw.triple.venn(
      area1 = area1,
      area2 = area2,
      area3 = area3,
      n12 = n12,
      n23 = n23,
      n13 = n13,
      n123 = n123,
      category = c("module",
                   "external",
                   "de"),
      lty = "blank",
      fill = c("red", "green", "blue")
    )
  }
  
  # Markers
  marker_names = intersect(module_names, de_names)
  
  if (!is.null(external_markers)) {
    marker_names = intersect(external_markers, marker_names)
  }
  
  return(sort(marker_names))
}