#' Function to caluclate cross-talk
#' @description `cross_talk()` calculates the cross-talk score for a provided matrix and vector of weights
#' @details The function calculates the cross-talk values for a provided matrix and a list of vectors of weights (fist row, second column).
#' The weights should be named by gene names
#' @param mat the matrix on which calculate the cross-talk. It should be the subset of the adjacency matrix of
#' a gene network, with on the rows the genes of a gene set, and on the columns the genes of anoher one to
#' @param weight the vector of weights associated to the gene in `mat`. The weights shoud be named by the genes
#' @return The function returns an array with:
#' \itemize{
#' \item ct: cross-talk value
#' \item ngenes_1, ngenes_2: number of gene in geneset 1 (rows) and geneset 2 (columns)
#' \item nlink: number of link between geneset 1 and 2
#' \item gene_1, gene_2: genes involved in the links between geneset 1 and 2
#' }
#' @importFrom collapse ss


cross_talk.perm <- function(mat=NULL, w1=NULL, w2=NULL){
  
  w1 <- ss(w1, rownames(mat))
  w2 <- ss(w2, colnames(mat))
  
  matW <- as.numeric(t(w1) %*% mat %*% w2)
  
  return(matW)
}
