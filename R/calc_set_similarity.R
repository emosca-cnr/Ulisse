#' Similarity between two sets
#' @param x,y two sets
#' @param method one among "jaccard", "overlap" and "average"
#' @return similarity



calc_set_similarity <- function(x, y, method=c("jaccard", "overlap", "average")){

  method <- match.arg(method)

  if(method=='jaccard' | method=='average' ){
    ssj <- sum(x %in% y)
    ssj <- ssj / length(unique(c(x, y)))
  }
  if(method=='overlap' | method=='average' ){
    sso <- sum(x %in% y)
    sso <- sso / min(length(x), length(y))
  }

  if(method=='jaccard')
    return(ssj)

  if(method=='overlap')
    return(sso)

  if(method=='average')
    return(mean(c(ssj, sso)))

}
