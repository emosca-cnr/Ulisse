#' Internal function

get_genes <- function(gene_set=NULL, wb=NULL){
  
  ans <- gene_set[gene_set %in% wb]
  if(!is.null(eg2sym)){
    ans <- sort(eg2sym$symbol[eg2sym$gene_id %in% ans])
  }
  
  ans <- paste0(ans, collapse = ";")
  return(ans)
  
}
