#' Function to compare results by using Spearman correlation
#' @description The function calculates the Spearman correlation between each pair of the provided results
#' @details The function takes as an input a list of elements, each  composed by a data.frame with three columns:
#' name of the first element (pathway or cell cluster), name of the second and the score that will be used for the comparison.
#' The function calculates the correlation score between each possible pairs and returns as an output a correlation matrix
#' where each rows and columns are named after the names of the elements of the list.
#' @param res_list a named list of data.frame, each composed by three columns: first and second elements of the cross-talk, 
#' score column.
#' @return The function returns a correlation matrix named after the names of `res_list`
#' @export



result_correlation <- function(res_list) {
  res_list <- lapply(res_list, function(x) {
    x[, 1:2] <- t(apply(x[, 1:2], 1, sort))
    colnames(x) <- c("ptw1", "ptw2", "score")
    return(x)
  })
  n <- names(res_list)
  res_list <- lapply(1:length(res_list), function(x) {
    tmp <- res_list[[x]]
    colnames(tmp)[3] <- paste0("score_", names(res_list)[x])
    return(tmp)
  })
  names(res_list) <- n
  
  res_mer <- Reduce(function(x, y) {
    merge(x, y, by =  c("ptw1", "ptw2"), all = T)
  }, res_list)
  # res_mer[,3:ncol(res_mer)] <- t(apply(res_mer[,3:ncol(res_mer)], 1, function(x) {
  #   x[is.na(x)] <- 0
  #   return(x)
  # }))
  cor.tab <- cor(res_mer[,3:ncol(res_mer)], method = "spearman", use = "na.or.complete")
  return(cor.tab)
}
