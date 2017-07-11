#' Reproducible Module analysis: Molecule Adjacency Matrix
#' 
#' Generate molecule-molecule adjacency matrix for molecule pairs exhibiting 
#' reproducible correlations across patients.
#' 
#' @param x an n by m numeric matrix, where m = num of samples, n = num of 
#'   molecules. Numerical values representing the expression level of the 
#'   molecule.
#' @param r factor. Sample replicate identity.
#' @param mad_constant scale factor passed to \code{\link[stats]{mad}} in 
#'   scaling z-transformed correlation coefficients.
#' @param var_thresh numeric. Molecules are only analyzed if their variance is 
#'   greater than this threshold in all replicate groups.
#' @param pair_pthresh numeric. Molecule pairs are called "reproducibly 
#'   correlated" if, across all replicate groups, their z-value corresponds to 
#'   a two-sided p-value below this threshold.
#'   
#' @export
#' 
#' @return Sparse adjacency matrix (\code{\link[Matrix]{dgCMatrix-class}}) with
#'   logical is_variable attribute flagging genes included in analysis.
#'   
#' @importFrom methods as
#' @importFrom stats var median mad pnorm cor
#'   
#' @examples
#' library("SummarizedExperiment")
#' data("fluidigm",package = "scRNAseq")
#' 
#' x = assay(fluidigm)[sample(x = seq_len(nrow(fluidigm)),size = 1000),]
#' r = factor(fluidigm$Coverage_Type)
#' A = get.repro.thresh.adjacency(x,r)
#' 
#' library("igraph")
#' g = graph.adjacency(A, mode = "undirected")
#' 
#' hist(degree(g))

get.repro.thresh.adjacency = function(x,
                                     r,
                                     mad_constant = 1.4826,
                                     var_thresh = 10 ^ -5,
                                     pair_pthresh = 0.01) {
  # Only test  that vary in all replicate pops
  is_var = rowSums(sapply(levels(r), function(p) {
    apply(x[, r == p], 1, stats::var)
  }) > var_thresh) == nlevels(r)
  
  vx = x[is_var, ]
  
  # All molecule pairs in all replicates
  pair_pmat = sapply(levels(r), function(p) {
    cor_mat = stats::cor(t(vx[, r == p]))
    pair_cor = atanh(cor_mat[upper.tri(cor_mat)])
    z_cor = (pair_cor - stats::median(pair_cor))
    z_cor = z_cor / stats::mad(pair_cor, constant = mad_constant)
    2 * stats::pnorm(-abs(z_cor))
  })
  
  # Matrix of reproducible pairs
  is_rep_pair = apply(pair_pmat < pair_pthresh, 1, all)
  rep_mat = matrix(
    0,
    nrow = nrow(vx),
    ncol = nrow(vx),
    dimnames = list(rownames(vx), rownames(vx))
  )
  rep_mat[upper.tri(rep_mat)] = is_rep_pair
  rep_mat = t(rep_mat)
  rep_mat[upper.tri(rep_mat)] = is_rep_pair
  
  rep_mat_full = matrix(
    0,
    nrow = nrow(x),
    ncol = nrow(x),
    dimnames = list(rownames(x), rownames(x))
  )
  rep_mat_full[rownames(vx), rownames(vx)] = rep_mat
  
  out_obj = methods::as(rep_mat_full, "dgCMatrix")
  
  attributes(out_obj)$is_variable = is_var
  
  out_obj
}

#' Reproducible Module analysis: Zero-Inflated Poisson Degree Distribution
#' 
#' Models degree distribution of an undirected graph by a zero-inflated Poisson
#' distribution, using the fitted Poisson parameters to define the null for
#' upper-tail probabilities.
#' 
#' @param g an undirected igraph graph (see
#'   \code{\link[igraph]{graph.adjacency}}).
#'   
#' @export
#' 
#' @return For each node, the upper tail probability under the fitted Poisson
#'   component distribution.
#'   
#' @importFrom igraph degree
#' @importFrom pscl zeroinfl
#' @importFrom stats ppois
#'   
#' @examples
#' library("SummarizedExperiment")
#' data("fluidigm",package = "scRNAseq")
#' 
#' x = assay(fluidigm)[sample(x = seq_len(nrow(fluidigm)),size = 1000),]
#' r = factor(fluidigm$Coverage_Type)
#' A = get.repro.thresh.adjacency(x,r)
#' 
#' library("igraph")
#' g = graph.adjacency(A, mode = "undirected")
#' 
#' pdegree = pzipdegree(g)
#' plot(degree(g),-log10(p.adjust(pdegree,method = "bonferroni")))

pzipdegree = function(g){
  d = igraph::degree(g)
  if(min(d) > 0){stop("invalid graph, minimum degree is not zero")}
  model.zip = pscl::zeroinfl(d ~ 1 | 1, dist = c("poisson"))
  stats::ppois(
    d,
    lambda = exp(model.zip$coefficients$count),
    lower.tail = FALSE
  )
}
  