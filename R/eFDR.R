#' empirical False Discovery Rate
#' @param real_values vactor of real (observed) values
#' @param all_values vactor of all values (real + permuted) values
#' @param mc.cores number of cores
#' @import parallel
#' @export
eFDR <- function(real_values, all_values, mc.cores=1){
  
  #FDR = (# x > permuted values / # x > real values) * (#real values / #permuted values)
  #approch used in GSEA
  
  fdr_values <- length(real_values) / length(all_values)
  fdr_values <- fdr_values*unlist(parallel::mclapply(real_values, function(x) sum(all_values >= x) / sum(real_values >= x), mc.cores = mc.cores))

  return(fdr_values)
  
}