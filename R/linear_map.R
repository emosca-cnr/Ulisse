#' Linear map
#'
#' @param x numeric vector
#' @param y.min min value of the new distribution
#' @param y.max max value of the new distribution
#' @export

linear_map <- function (x, y.min, y.max){
  m <- ((y.max - y.min)/(max(x) - min(x)))
  c <- y.max - m * max(x)
  y <- m * x + c
  return(y)
}
