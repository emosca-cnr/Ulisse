#' Caluclate communication between two clusters
#' @param DEG_list a gene list composed by a vector for both clusters. The vector should be a weight vector named by gene
#' @param gene_network_adj gene network adjacency matrix
#' @param ligand vector of ligand genes
#' @param receptor vector of receptor genes
#' @param mc_cores_perm number of threads to be used to calcualte permutation
#' @param k number of permutation of the adjacency matrix
#' @import parallel
#' @export


cluster_communication <- function(DEG_list, gene_network_adj, k = 9, ligand = ligand, receptor = receptor, 
                                  mc_cores_perm = 1) {
  ligand <- ligand[which(ligand %in% rownames(gene_network_adj))]
  receptor <- receptor[which(receptor %in% rownames(gene_network_adj))]
  mixed <- intersect(ligand, receptor)
  if(length(mixed) > 0) {
    ligand <- ligand[-which(ligand %in% mixed)]
    receptor <- receptor[-which(receptor %in% mixed)]
  }
  
  gene_network_adj <- sign(gene_network_adj)
  gene_network_adj[ligand, ligand] <- 0
  gene_network_adj[receptor, receptor] <- 0
  
  perm_list <- mclapply(1:k, function(x) {
    tmp <- matrix(as.numeric(gene_network_adj), ncol = ncol(gene_network_adj), 
                  dimnames = list(sample(rownames(gene_network_adj), nrow(gene_network_adj))))
    colnames(tmp) <- rownames(tmp)
    return(tmp)
  }, mc.cores = mc_cores_perm
  )
  perm_list <- c(list(gene_network_adj), perm_list)
  ans <- mclapply(1:(k+1), function(x) {
    mat <- perm_list[[x]][names(DEG_list[[1]]), names(DEG_list[[2]]), drop = F]
    matW <- t(DEG_list[[1]]) %*% as.matrix(mat)
    matW <- matW %*% (DEG_list[[2]])
    matW <- data.frame(cl1 = names(DEG_list)[1],
                       cl2 = names(DEG_list)[2],
                       cct = matW,
                       ngenes_cl1 = length(DEG_list[[1]]),
                       ngenes_cl2 = length(DEG_list[[2]]),
                       nlink = sum(mat),
                       stringsAsFactors = F)
    
    return(matW)
  }, mc.cores = mc_cores_perm)
  
  
  p_list <- mclapply(ans, function(x) {
    rownames(x) <- paste(x$cl1, x$cl2, sep = "|")
    x <- x[, "cct", drop = F]
    return(x)
  }, mc.cores = mc_cores_perm)
  p_val <- calc_p(p_list)
  
  out <- ans[[1]]
  out$p_value <- p_val
  out$eFDR <- eFDR(real_values = as.vector(unlist(out$cct)), 
                   all_values = as.vector(unlist(p_list)))
  
  mat <- gene_network_adj[names(DEG_list[[1]]), names(DEG_list[[2]]), drop=F]
  idx <- which(mat==1, arr.ind = T)
  cl_ct <- data.frame(cl1 = rep(names(DEG_list)[1], nrow(idx)),
                      cl1_gene = rownames(idx),
                      cl1_type = rep("ligand", nrow(idx)),
                      cl2 = rep(names(DEG_list)[2], nrow(idx)),
                      cl2_gene = colnames(mat)[idx[,2]],
                      cl2_type = rep("ligand", nrow(idx)),
                      stringsAsFactors = F)
  cl_ct$cl1_type[which(cl_ct$cl1_gene %in% receptor)] <- "receptor"
  cl_ct$cl2_type[which(cl_ct$cl2_gene %in% receptor)] <- "receptor"
  cl_ct$cl1_type[which(cl_ct$cl1_gene %in% mixed)] <- "ligand/receptor"
  cl_ct$cl2_type[which(cl_ct$cl2_gene %in% mixed)] <- "ligand/receptor"
  
  
  return(list(communications_info = cl_ct, cc_communications= out))
}
