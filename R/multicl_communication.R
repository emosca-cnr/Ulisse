#' Caluclate communication between multiple clusters
#' @param DEG_list a list composed by gene list for each couple of clusters. Each gene list should be composed by weight vector named by gene for each cluster
#' @param gene_network_adj gene network adjacency matrix
#' @param ligand vector of ligand genes
#' @param receptor vector of receptor genes
#' @param mc_cores_perm number of threads to be used to calcualte permutation
#' @param k number of permutation of the adjacency matrix
#' @import parallel
#' @export

multicl_communication <- function(DEG_list, gene_network_adj, k = 9, ligand = ligand, receptor = receptor, 
                                  mc_cores_perm = 1) {
  perm_list <- mclapply(1:k, function(x) {
    tmp <- matrix(as.numeric(gene_network_adj), ncol = ncol(gene_network_adj), 
                  dimnames = list(sample(rownames(gene_network_adj), nrow(gene_network_adj))))
    colnames(tmp) <- rownames(tmp)
    return(tmp)
  }, mc.cores = mc_cores_perm
  )
  perm_list <- c(list(gene_network_adj), perm_list)
  ans <- mclapply(1:(k+1), function(x) {
    cl_ct <- lapply(1:length(DEG_list), function(i) {
      matW <- cluster_communication(DEG_list = DEG_list[[i]], gene_network_adj = perm_list[[x]], 
                                    ligand = ligand, receptor = receptor, k=0, mc_cores_perm = mc_cores_perm)
      return(matW)
    })
    
    
    return(cl_ct)
  }, mc.cores = mc_cores_perm)
  
  comm_list <- list()
  for(j in 1:length(ans)) {
    tmp <- ans[[j]]
    tmp <- lapply(tmp, function(x) return(x[[2]]))
    comm_list[[j]] <- do.call(rbind, tmp)
  }
  comm_info <- do.call(rbind, lapply(ans[[1]], function(x) return(x[[1]])))
  
  p_list <- mclapply(comm_list, function(x) {
    rownames(x) <- paste(x$cl1, x$cl2, sep = "|")
    x <- x[, "cct", drop = F]
    return(x)
  }, mc.cores = mc_cores_perm)
  p_val <- calc_p(p_list)
  
  out <- comm_list[[1]]
  out$p_value <- p_val
  out$eFDR <- eFDR(real_values = as.vector(unlist(out$cct)), 
                   all_values = as.vector(unlist(p_list)))
  return(list(communications_info = comm_info, cc_communications= out))
}



