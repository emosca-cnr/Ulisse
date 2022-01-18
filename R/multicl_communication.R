#' Function to caluclate communication between multiple clusters
#' @description multicl_communication calculates a score for the communication between the cluster pairs provided
#' @details The `gene_network_adj` is modified by using `ligand` and `receptor` information to remove the  
#' ligand-ligand and receptor-receptor interactions. Then, CT formula is used to calculate a score for the communication 
#' of each cluster pair considering only the genes listed in the `DEG_list` and their weights. The same approach is 
#' applied also to the `k` permuted version of `gene_network_adj` to calculate the p-value and the FDR
#' @param DEG_list a gene list composed by a list for each cluster pair, composed by vector for both clusters. 
#' The vector should be a weight vector named by gene
#' @param gene_network_adj gene network adjacency matrix
#' @param ligand vector of ligand genes
#' @param receptor vector of receptor genes
#' @param mc_cores_perm number of threads to be used to calcualte permutation
#' @param k number of permutation of the adjacency matrix
#' @return The function returns a list of two objects:
#' \enumerate{
#' \item communications_info: a data.frame with the detailed description of the communication between the clusters
#' \itemize{
#' \item cl1, cl2: name of the clusters (derived from the name of the list in `DEG_list`)
#' \item cl1_gene, cl2_gene: genes that are ligand or receptor in cl1 or cl2 that partecipate to the communication
#' \item cl1_type, cl2_type: type of the gene in cl1/cl2 that participate in the communication. Can be ligand, receptor
#' or ligand/receptor if the gene is present in both `ligand` and `receptor` vectors and it is supposed it behave 
#' as a ligand 
#' in some conditions and as a receptors in others
#' \item score: the score of the communication between the two genes in cl1 and cl2, claculated considering their weights
#' }
#' \item {cc_communications}: a data.frame with the communication score
#' \itemize{
#' \item cl1, cl2: names of the clusters
#' \item ccc_score: communication score calculated on all the interacting ligand-receptor genes in the cluster pair
#' \item ngenes_cl1, ngenes_cl2: number of genes in cl1 and cl2, respectively, that participate to the communication score
#' \item nlink: number of links beteewn the genes in cl1 and cl2
#' \item p_value: empirical p-value calculated by using the permutation approcah
#' \item eFDR: empirical FDR
#' }
#' }
#' @examples 
#' cl1 <- rep(1, 4)
#' names(cl1) <- c( "A", "B", "C")
#' cl2 <- rep(1, 4)
#' names(cl2) <- c("D", "E", "F")
#' cl3 <- rep(1, 4)
#' names(cl3) <- c("G", "H", "I")
#' DEG_list <- list("cl1|cl2" = list("cl1" = cl1, "cl2" = cl2), 
#' "cl2|cl3" = list("cl2" = cl2, "cl3" = cl3),
#' "cl1|cl3" = list("cl1" = cl1, "cl3" = cl3))
#' adj <- matrix(data = sample(c(0,1), 9*9, replace = T), nrow = 9, ncol = 9, 
#' dimnames = list(LETTERS[1:9], LETTERS[1:9]))
#' ligand <- c("A", "D", "G", "F", "C", "I")
#' receptor <- c("B", "E", "H", "C", "I")
#' cc <- multicl_communication(DEG_list = DEG_list, gene_network_adj = adj, ligand = ligand, 
#' receptor = receptor)
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
    x <- x[, "ccc_score", drop = F]
    return(x)
  }, mc.cores = mc_cores_perm)
  p_val <- calc_p(p_list)
  
  out <- comm_list[[1]]
  out$p_value <- p_val
  out$eFDR <- eFDR(real_values = as.vector(unlist(out$ccc_score)), 
                   all_values = as.vector(unlist(p_list)))
  return(list(communications_info = comm_info, cc_communications= out))
}




