#' Function to caluclate communication between cell clusters
#' @description cluster_communication calculates a score for the communication between all the possible pairs 
#' in the cluster list provided
#' @details Firstly, `k` permuted version of  `gene_network_adj` are created. Then, CT formula is used to 
#' calculate a score for the communication between all possible cluster pairs. The same approach is 
#' applied also to the `k` permuted version of `gene_network_adj` to calculate the p-value and FDR
#' @param cl_list the output of `preparing_cl_list()` or a gene list composed by a vector for each cluster where each vector is
#' a named gene rank or weight vector.
#' @param gene_network_adj gene network adjacency matrix
#' @param mc_cores_perm number of threads to be used to calcualte permutation
#' @param mc_cores_ccc number of threads to be used for CCC calculation
#' @param k number of permutation of the adjacency matrix
#' @return The function returns a list of two objects:
#' \enumerate{
#' \item communications_info: a data.frame with the detailed description of the communication between the two clusters
#' \itemize{
#' \item cl1, cl2: name of the two clusters (derived from the name of the list in `DEG_list`)
#' \item cl1_gene, cl2_gene: genes that are ligand or receptor in cl1 or cl2 that partecipate to the communication
#' \item score: the score of the communication between the two genes in cl1 and cl2, claculated considering their weights
#' }
#' \item {cc_communications}: a data.frame with the communication score
#' \itemize{
#' \item cl1, cl2: names of the two clusters
#' \item ccc_score: communication score calculated on all the interacting ligand/receptor genes in the two clusters
#' \item ngenes_cl1, ngenes_cl2: number of genes in cl1 and cl2, respectively, that participate to the communication score
#' \item nlink: number of links beteewn the genes in cl1 and cl2
#' \item p_value: empirical p-value calculated by using the permutation approcah
#' \item eFDR: false discovery rate (FDR)
#' \item genes_cl1, genes_cl2: list of the genes involved in the communication in the clusters
#' }
#' }
#' @examples 
#' \dontrun{
#' #' cl1 <- rep(1, 4)
#' names(cl1) <- c( "A", "B", "C", "D")
#' cl2 <- rep(1, 4)
#' names(cl2) <- c( "E", "F", "G", "H")
#' cl_list <- list("cl1" = cl1, "cl2" = cl2)
#' adj <- matrix(data = sample(c(0,1), 8*8, replace = TRUE), nrow = 8, ncol = 8, 
#' dimnames = list(LETTERS[1:8], LETTERS[1:8]))
#' cc <- cluster_communication(cl_list = cl_list, gene_network_adj = adj)
#' }
#' @import igraph
#' @import parallel
#' @import stringi
#' @export

cluster_communication <- function(cl_list, gene_network_adj=g.adj, k = 9, 
                                  mc_cores_perm = 1, mc_cores_ccc = 1) {
  
  gene_network_adj <- sign(gene_network_adj)
  
  perm_list <- parallel::mclapply(1:k, function(x) {
    tmp <- gene_network_adj
    rownames(tmp) <- sample(rownames(tmp), nrow(tmp))
    colnames(tmp) <- rownames(tmp)
    return(tmp)
  }, mc.cores = mc_cores_perm
  )
  perm_list <- c(list(gene_network_adj), perm_list)
  
  comb_p <- expand.grid(names(cl_list), names(cl_list))
  g <- graph_from_edgelist(as.matrix(comb_p), directed = F)
  g <- igraph::simplify(g, remove.multiple = T, remove.loops = T)
  comb_p <- as_edgelist(g)
  
  ans <- parallel::mclapply(1:(k+1), function(x) {
    tmp <- parallel::mclapply(1:nrow(comb_p), function(z) {
      tab <- perm_list[[x]]
      cl1 <- cl_list[[comb_p[z, 1]]]
      ncl1 <- names(cl1)
      cl2 <- cl_list[[comb_p[z, 2]]]
      ncl2 <- names(cl2)
      mat <- tab[ncl1, ncl2, drop = F]
      matW <- t(cl1) %*% as.matrix(mat)
      matW <- matW %*% (cl2)
      
      row.col.idx <- which(mat == 1, arr.ind = T)
      row.n <- rownames(mat)[row.col.idx[,1]]
      row.n <- ncl1[ncl1 %in% row.n]
      col.n <- colnames(mat)[row.col.idx[,2]]
      col.n <- ncl2[ncl2 %in% col.n]
      
      tab.out <- array(data = NA, dim=10)
      tab.out[1] <- comb_p[z, 1]
      tab.out[2] <- comb_p[z, 2]
      tab.out[3] <- matW
      tab.out[4] <- length(row.n)
      tab.out[5] <- length(col.n)
      tab.out[6] <- sum(mat)
      tab.out[7] <- stri_c(row.n, collapse = ";")
      tab.out[8] <- stri_c(col.n, collapse = ";")
      tab.out[9] <- sum(cl1)
      tab.out[10] <- sum(cl2)
      
      
      
      return(tab.out)
    },mc.cores = mc_cores_ccc) 
    tmp <- do.call(rbind, tmp)
    colnames(tmp) <- c("cl1", "cl2", "ccc_score", "ngenes_cl1", "ngenes_cl2",
                       "nlink", "genes_cl1", "genes_cl2", "gene_weight_cl1", "gene_weight_cl2")
    return(tmp)
  }, mc.cores = mc_cores_perm)
  
  p_list <- parallel::mclapply(ans, function(x) {
    rownames(x) <- paste(x[,1], x[,2], sep = "|")
    x <- x[, "ccc_score", drop = F]
    return(x)
  }, mc.cores = mc_cores_perm)
  p_val <- calc_p(p_list)
  
  out <- ans[[1]]
  eFDR <- Ulisse::eFDR(real_values = as.vector(unlist(p_list[[1]])), all_values = as.vector(unlist(p_list)))
  out <- data.frame(out, stringsAsFactors = F)
  out$p_value <- as.vector(p_val)
  out$FDR <- eFDR
  
  ct_info <- parallel::mclapply(1:nrow(comb_p), function(z) {
    cl1 <- cl_list[[comb_p[z, 1]]]
    ncl1 <- names(cl1)
    cl2 <- cl_list[[comb_p[z, 2]]]
    ncl2 <- names(cl2)
    mat <- gene_network_adj[ncl1, ncl2, drop = F]
    
    row.col.idx <- which(mat == 1, arr.ind = T)
    
    score = cl1[rownames(row.col.idx)] * cl2[colnames(mat)[row.col.idx[,2]]]
    cl1 = rep(comb_p[z, 1], nrow(row.col.idx))
    cl1_gene = rownames(row.col.idx)
    cl2 = rep(comb_p[z, 2], nrow(row.col.idx))
    cl2_gene = colnames(mat)[row.col.idx[,2]]
    
    cl_ct <- matrix(c(cl1, cl1_gene, cl2, cl2_gene, score), nrow = nrow(row.col.idx))
    return(cl_ct)
  },mc.cores = mc_cores_ccc)
  ct_info <- do.call(rbind, ct_info)
  colnames(ct_info) <- c("cl1", "cl1_gene", "cl2", "cl2_gene", "score")
  ct_info <- data.frame(ct_info, stringsAsFactors = F)
  
  return(list(communications_info = ct_info, cc_communications= out))
}



