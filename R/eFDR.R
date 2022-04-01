#' empirical False Discovery Rate
#' @param real_values vactor of real (observed) values
#' @param all_values vactor of all values (real + permuted) values
#' @param mc.cores number of cores
#' @param correct.max logical, to decide if to correct the maximum values to 1 or not
#' @import parallel
#' @export
eFDR <- function(real_values, all_values, mc.cores=1, correct.max = TRUE){
  
  #FDR = (# x > permuted values / # x > real values) * (#real values / #permuted values)
  #approch used in GSEA
  
  fdr_values <- length(real_values) / length(all_values)
  fdr_values <- fdr_values*unlist(parallel::mclapply(real_values, function(x) sum(all_values >= x) / sum(real_values >= x), mc.cores = mc.cores))
  if(correct.max) {
    fdr_values[fdr_values>1] <- 1
  }

  return(fdr_values)
  
}