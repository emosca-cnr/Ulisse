#' Returns filtered msigdb geneset in a format compatible with cross-talks calculations
#' @description Download MsidgDB data, filter them and prepare them in the format needed by cross-talks calculation
#' @details The function downloads the pathway data for the `species` of interest. The selection of the database to 
#' download is made by using the parameter `category` and `subcategory`. The data downloaded are then filtered to
#' maintain only the `genes` of interest and to remove pathway with a number of genes lower of `min_size` and higher of
#' `max_size`
#' @param species Species name
#' @param category MSigDB category name, see `msigdbr_collections()`
#' @param subcategory MsigDB subcategory name, see `msigdbr_collections()`
#' @param type Gene name of interest, can be `gene_symbol` or `entrez_gene`
#' @param genes vector with gene gene name of interest
#' @param min_size minimum number of gene present in each gene set to be considered
#' @param max_size maximum number of gene present in each gene set to be considered
#' @return The output is a list of two output:
#' \enumerate{
#' \item msigdb_output: a data.frame with the data downloaded by MsigDB
#' \item path_list: a list of pathway. Each list is composed by a vector with the names of the genes of interest 
#' that are part of that pathway
#' }
#' @import msigdbr
#' @export

pathway_data <- function (species, category = NULL, subcategory = NULL, type = "gene_symbol", genes,
                          min_size=1, max_size=500 ) {
  msig_out <- msigdbr::msigdbr(species = species, category = category, 
                           subcategory = subcategory)
  if (type == "gene_symbol") {
    msig <- msig_out[, c("gene_symbol", "gs_name")]
    msig <- unique(msig)
    msig_list <- split(x = msig$gene_symbol, f = msig$gs_name)
  } else if (type == "entrez_gene") {
    msig <- msig_out[, c("entrez_gene", "gs_name")]
    msig <- unique(msig)
    msig_list <- split(x = msig$entrez_gene, f = msig$gs_name)
  }
  output <- list(msigdb_output = msig_out, path_list = msig_list)
  
  output$path_list <- Ulisse::filter_gsl(gsl = output$path_list, universe = genes, 
                                         min_size = min_size, max_size = max_size)
  
  return(output)
}
