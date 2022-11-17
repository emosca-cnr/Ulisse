#' Function to caluclate cross-talk
#' @description `cross_talk()` calculates the cross-talk score for a provided matrix and vector of weights
#' @details The function calculates the cross-talk values for a provided matrix and a vector of weights.
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
#' @importFrom stringi stri_c
#' @export


cross_talk <- function(mat, weight) {
  g.1 <- as.character(rownames(mat))
  g.2 <- as.character(colnames(mat))
  wg.1 <- weight[["g1"]]
  wg.2 <- weight[["g2"]]
  matW <- t(wg.1[g.1]) %*% mat
  matW <- matW %*% wg.2[g.2, drop = F]
  row.col.idx <- which(mat == 1, arr.ind = T)
  row.n <- funique(g.1[row.col.idx[,1]])
  row.n <- row.n[row.n %in% names(wg.1)]
  col.n <- funique(g.2[row.col.idx[,2]])
  col.n <- col.n[col.n %in% names(wg.2)]
  
  #mat.out <- matrix(data = NA, nrow = 1, ncol = 8)
  mat.out <- array(data = NA, dim = 8, dimnames = list(c("ct", "ngenes_1",
                                                         "ngenes_2", "nlink","weight_1", 
                                                         "weight_2", "gene_1", "gene_2")))
  mat.out[1] = matW
  mat.out[2] = length(row.n)
  mat.out[3] = length(col.n)
  mat.out[4] = sum(mat)
  mat.out[5] = sum(wg.1)
  mat.out[6] = sum(wg.2)
  mat.out[7] = stri_c(row.n, collapse = ";")
  mat.out[8] = stri_c(col.n, collapse = ";")
  
  return(mat.out)
  
}