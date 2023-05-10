#' Estimation of p values
#' @description This function is used to calculate p-value from a list of matrices
#' @details the functions calculate p-value for a permutation-based approach. It takes as an input a list of matrices, where the first is the one with the real data
#' and the other are the data obtained through permutations. Then, the function searches for how many values in the permutations are equal or higher to the real.
#' These values are then divided by the number matrices (1 + number of permutations)
#' @param X list of matrices, where the first is the one obtained with real data and the others are obtained through permutations. 
#' The matrices should be 1 column matrices with the same order.
#' @export
#'
calc_p <- function(X){
  p <- matrix(0, nrow=nrow(X[[1]]), ncol=ncol(X[[1]]), 
              dimnames = list(rownames(X[[1]]), colnames(X[[1]])))
  for(i in 1:length(X)){
    idx <- X[[i]] >= X[[1]]
    p[idx] <- p[idx] + 1
  }
  p <- p / length(X)
  return(p)
}