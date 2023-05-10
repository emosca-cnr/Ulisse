#' Build list for cross-talk calculation from MSigDB database list
#' @description The function build a gene set list compatible with `gs_cross-talk()` for cross-talk calculation using a desired MSigDB (using msigdbr package)
#' @details `preparing_msigdb_list()` uses a MSigDB database to build the list needed by `gs_cross_talk()` for cross-talk calculation. In details, it calls
#'  `pathway_data()` function to download the database by using msigdbr package to then build the gene-set list
#' @param species Species name
#' @param category MSigDB category name, see `msigdbr_collections()`
#' @param subcategory MsigDB subcategory name, see `msigdbr_collections()`
#' @param type Gene name of interest, can be `gene_symbol` or `entrez_gene`
#' @param genes vector with gene gene name of interest
#' @param weights vector with `genes` weights in the same order. If `NULL` the function assigns `1` to each gene
#' @param min_size,max_size filtering the gene sets to maintain only the ones that have at least `min_size` or less than `max_size` number of genes
#' @return The output is a gene-set list with a vector for each gene-set, composed by weights named after the genes
#' @importFrom stats setNames
#' @export

preparing_msigdb_list <- function (species, category = NULL, subcategory = NULL, type = "gene_symbol", genes, weights,
                                   min_size=1, max_size=500) {
  if(is.null(weights)) {
    weights <- rep(1, length(genes))
  }
  weights <- setNames(weights, genes)
  ptw <- Ulisse::pathway_data(species = species, category = category, subcategory = subcategory, type = type, genes = genes,
                            min_size= min_size, max_size = max_size)
  db_out <- lapply(ptw$path_list, unique)
  db_out <- lapply(db_out, function(x) x <- weights[x])
  
  return(db_out)
}
