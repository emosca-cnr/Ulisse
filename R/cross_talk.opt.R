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
#' @importFrom stringi stri_c
#' @importFrom collapse ss
#' @import Matrix


cross_talk.opt <- function(mat=NULL, w1=NULL, w2=NULL) {
  
  w1 <- ss(w1, rownames(mat))
  w2 <- ss(w2, colnames(mat))

  nE <- sum(mat)
  
  matW <- as.numeric(t(w1) %*% mat %*% w2)

  real_link <- as.numeric(sign(t(w1)) %*% mat %*% sign(w2))
  
  #mat.sub <- mat[ names(w1)[!w1 == 0], names(w2)[!w2 == 0], drop = F]
  mat.sub <- ss(mat, rownames(w1)[w1 != 0], rownames(w2)[w2 != 0])
  row.col.idx <- which(mat.sub == 1, arr.ind = T)
  
  if(length(row.col.idx) > 0) {
    row.n <- funique(rownames(mat.sub)[row.col.idx[, 1]])
    #row.n <- row.n[row.n %in% names(wg.1)[!wg.1 == 0]]
    col.n <- funique(colnames(mat.sub)[row.col.idx[, 2]])
    #col.n <- col.n[col.n %in% names(wg.2)[!wg.2 == 0]]
  } else {
    row.n <- NULL
    col.n <- NULL
  }
  
  mat.out <- data.frame(
    c = matW,
    #r = matW/nE,
    S1_size = length(w1),
    S2_size = length(w2),
    S1_S2_size = length(row.n),
    S2_S1_size = length(col.n),
    dL = real_link,
    L = nE,
    r_c = real_link/nE,
    u1 = sum(w1[row.n, ]),
    u2 = sum(w2[col.n, ]),
    S1 = stri_c(row.n, collapse = ";"),
    S2 = stri_c(col.n, collapse = ";"),
    stringsAsFactors=FALSE
  )
  
  return(mat.out)
  
}