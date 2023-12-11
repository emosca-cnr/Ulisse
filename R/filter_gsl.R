#' Filter a gene set list
#' @param gsl gene set list (named list)
#' @param universe set of all possible values for items of gsl
#' @param min_size minimum gene set size
#' @param max_size maximum gene set size
#' @export


filter_gsl <- function(gsl=NULL, universe=NULL, min_size=5, max_size=500){
  
  ans <- gsl
  
  gs_size <- lengths(ans)
  idx_keep <- gs_size >= min_size & gs_size <= max_size
  ans <- ans[idx_keep]
  
  ans <- lapply(gsl, function(x) x[x %in% universe])
  ans <- lapply(ans, unique)
  
  gs_size <- lengths(ans)
  idx_keep <- gs_size >= min_size & gs_size <= max_size
  
  ans <- ans[idx_keep]
  
  return(ans)
  
}