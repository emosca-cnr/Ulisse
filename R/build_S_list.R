#' Build list for cross-talk calculation
#' @description The function build a gene-set list compatible with `gs_cross-talk()` for cross-talk calculation starting from a provided gene set database
#'  and genes (with weights) of interest.
#' @details `build_S_list()` uses a gene-set database data to build the list needed by `gs_cross_talk()` for cross-talk calculation.
#' @param S_tab a table composed by two (mandatory) or three columns: the first contain the names of the gene sets, the second the name of the genes. The third column
#'  (if provided) should be the gene weights. If `ncol(S_tab) == 3`, `g_u` param is ignored
#'  weights in the gene sets (i.e.: cell-cell communication), there should be a third column with the weights of each gene in the gene sets.
#' @param g_u named vector. Should contain the genes (as names) and their weights. Ignored if `ncol(S_tab) == 3`.
#' @param min_size,max_size filtering the gene sets to maintain only the ones that 
#'  have at least `min_size` or less than `max_size` number of genes. Set to `NULL` to skip filtering
#' @param universe genes in the adjacency matrix that will be used in `gs_cross_talk()`. Set to `NULL` to skip filtering
#' @return The output is a gene-set list with a vector for each gene-set, composed by weights named after the genes
#' @importFrom stats setNames
#' @export

build_S_list <- function (S_tab, g_u, universe, min_size=1, max_size=500) {
  if(ncol(S_tab) == 3) {
    out <- split(S_tab, S_tab[, 1])
    out <- lapply(out, function(tab) tab <- setNames(tab[, 3], tab[, 2]))
  } else {
    out <- split(S_tab[, 2], S_tab[, 1])
    out <- lapply(out, function(l) l <- unique(l))
    out <- lapply(out, function(l) l <- g_u[as.character(l)])
  }
  
  if(!is.null(universe)) {
    out <- lapply(out, function(l) l <- l[names(l) %in% universe])
  }
  
  if(!is.null(min_size) & !is.null(max_size)) {
    idx <- lengths(out) >= min_size & lengths(out) <= max_size
    out <- out[idx]
  } else if(!is.null(min_size) & is.null(max_size)) {
    idx <- lengths(out) >= min_size 
    out <- out[idx]
  } else if(is.null(min_size) & !is.null(max_size)) {
    idx <- lengths(out) <= max_size
    out <- out[idx]
  } 
  
  
  
  return(out)
}
