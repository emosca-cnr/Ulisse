#' Internal function of pathtools
#' @param rll list of named ranked vectors
#' @param perm vector of permuted names
#' @param gs gene set
calc_gs_perm <- function(rll=NULL, perm=NULL, gs=NULL){
  
  out <- unlist(lapply(rll, function(x) es(idx = which(perm %in% gs), array(x, dimnames=list(perm)))[, 1]))

  return(out)

}