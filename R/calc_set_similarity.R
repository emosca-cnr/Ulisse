#' The function calculates similarities between two sets.
#' @details `calc_set_similarity` calculates similarities between two sets by using `jaccard`, `overlap`, `average` methods.
#' @description The function takes as input two vertices and calculates their similarity by using three different methods:
#' \itemize{
#' \item Jaccard = number of common elements between x and y, divided by the number of unique elements present in c(x, y)
#' \item Overlap = number of common elements between x and y. divided by the minimum length of unique elements present in either x or y
#' \item Average = the mean between Jaccard and Overlap measure
#' }
#' @param x,y two sets provided as vectors
#' @param method one among "jaccard", "overlap" and "average"
#' @return similarity


calc_set_similarity <- function(x, y, method=c("jaccard", "overlap", "average"))
  {
  method <- match.arg(method)
  if(method == 'jaccard' | method == 'average' ){
    ssj <- sum(x %in% y)
    ssj <- ssj / length(unique(c(x, y)))
  }
  if(method == 'overlap' | method == 'average' ){
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