#' Function to prepare cluster list for `cluster_communication()`
#' @description `preparing_cl_list()` produce the cluster list needed for CCC calculation
#' @details The function is used to prepare the cluster list needed as an input in `cluster_communication()`. In detail,
#' the function needs the gene_by_cell matrix used in single cell analysis. We suggest to use the normalized data to avoid
#' accounting for differences in gene counts. The matrix is binarized by assigning a `1` to all the genes that have a normalized 
#' expression value equal or higher than `mean_t`, `0` otherwise. Then, for each gene in cluster, the function discards the 
#' gene expressed in a number of cell lower then `cell_t` and cacluates the mean frequency of over-threshold 
#' expression between all cells of the cluster. The function returns a vector for each cluster, each one composed by 
#' the frequencies of the genes and named after the respecitve genes. The vectors are filtered to maintain only the
#' genes in `universe`, that are the gene of interest
#' @param mtx the gene-by-cell matrix used in single cell data analysis
#' @param clusters an  vector with cluster membership of all the cells in `mtx`
#' @param universe names of the genes of interest. Used to filter the cluster gene list
#' @param mean_t threshold on the expression used to binarize the mexpression matrix. This value should be decided considering 
#' the distribution of the mean expression of all the genes (without considering the zeros)
#' @param cell_t minumum number of cells in which a gene should be expressed over `mean_t` threshold in a cluster
#' to be considered in subsequent analysis
#' @return The function returns a gene list composed by a vector for each cluster provided. Each vector is composed by the
#' frequency of over-threshold expression of the genes in the cluster, named by the respective genes
#' @export


preparing_cl_list <- function(mtx, clusters, mean_t=1, 
                              cell_t = 5, universe) {
  mtx <- as.matrix(mtx)
  
  mtx[mtx < mean_t] <- 0
  mtx[mtx >= mean_t] <- 1
  
  cl <- unique(clusters)
  cl_list <- list()
  for (i in cl) {
    cell <- colnames(mtx)[clusters == i]
    tmp <- mtx[, cell]
    tmp <- tmp[rowSums(tmp) >= cell_t,]
    tmp2 <- rowMeans(tmp)
    names(tmp2) <- rownames(tmp)
    cl_list <- c(cl_list, list(tmp2))
  }
  names(cl_list) <- cl
  cl_list <- lapply(cl_list, function(x) x[names(x) %in% universe])
  return(cl_list)
}

