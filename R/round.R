#' Scientific approximation
#'
#' @param x numeric
#' @param y number of digits
#' @export

round <- function(x, y=0){
  return(trunc(x*10^y + 0.5)/(10^y))
}
