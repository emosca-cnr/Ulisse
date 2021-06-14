#' Estimation of p values
#' @param X list of matrices, where the first is the one obtained with real data
#' @export
#'
calc_p <- function(X){

  p <- matrix(0, nrow=nrow(X[[1]]), ncol=ncol(X[[1]]), dimnames = list(rownames(X[[1]]), colnames(X[[1]])))

  for(i in 1:length(X)){
    idx <- X[[i]] >= X[[1]]
    p[idx] <- p[idx] + 1
  }

  p <- p / length(X)

  return(p)


}
