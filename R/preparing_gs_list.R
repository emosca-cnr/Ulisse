#' Build list for cross-talk calculation from database list
#' @description The function build a gene-set list compatible with `gs_cross-talk()` for cross-talk calculation starting from a provided gene set database
#'  and genes (with weights) of interest.
#' @details `preparing_gs_list()` uses a gene-set database data to build the list needed by `gs_cross_talk()` for cross-talk calculation.
#' @param gs_names names of the gene-sets, to be used to build the gene-set list. Should have the same length of `gs_gene` 
#' @param gs_genes names of the genes in the gene-sets `gs_genes`
#' @param weights named vector with a weight for each gene in `gs_gene`. If `NULL` the function assigns `1` to each gene
#' @return The output is a gene-set list with a vector for each gene-set, composed by weights named after the genes
#' @importFrom stats setNames
#' @export

preparing_gs_list <- function (gs_names , gs_genes, weights, min_size=1, max_size=500) {
  if(is.null(weights)) {
    genes <- unique(gs_genes)
    weights <- setNames(rep(1, length(genes)), genes)
  }
  
  db_list <- split(gs_genes, gs_names)
  db_list <- lapply(db_list, unique)
  db_list <- lapply(db_list, function(x) x <- x[x %in% names(weights)])
  db_list <- db_list[lengths(db_list)>0]
  db_list <- Ulisse::filter_gsl(gsl = db_list, 
                                         universe = gs_genes, min_size = min_size, max_size = max_size)
  db_list <- lapply(db_list, function(x) x <- weights[x])
  return(db_list)
}
