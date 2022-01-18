#' Function to caluclate communication between two clusters
#' @description cluster_communication calculates a score for the communication between the two provided clusters
#' @details The `gene_network_adj` is modified by using `ligand` and `receptor` information to remove the  
#' ligand-ligand and receptor-receptor interactions. Then, CT formula is used to calculate a score for the communication 
#' between the two clusters considering only the genes listed in the `DEG_list` and their weights. The same approach is 
#' applied also to the `k` permuted version of `gene_network_adj` to calculate the p-value.
#' @param DEG_list a gene list composed by a vector for both clusters. 
#' The vector should be a weight vector named by gene
#' @param gene_network_adj gene network adjacency matrix
#' @param ligand vector of ligand genes
#' @param receptor vector of receptor genes
#' @param mc_cores_perm number of threads to be used to calcualte permutation
#' @param k number of permutation of the adjacency matrix
#' @return The function returns a list of two objects:
#' \enumerate{
#' \item communications_info: a data.frame with the detailed description of the communication between the two clusters
#' \itemize{
#' \item cl1, cl2: name of the two clusters (derived from the name of the list in `DEG_list`)
#' \item cl1_gene, cl2_gene: genes that are ligand or receptor in cl1 or cl2 that partecipate to the communication
#' \item cl1_type, cl2_type: type of the gene in cl1/cl2 that participate in the communication. Can be ligand, receptor
#' or ligand/receptor if the gene is present in both `ligand` and `receptor` vectors and it is supposed it behave 
#' as a ligand 
#' in some conditions and as a receptors in others
#' \item score: the score of the communication between the two genes in cl1 and cl2, claculated considering their weights
#' }
#' \item {cc_communications}: a data.frame with the communication score
#' \itemize{
#' \item cl1, cl2: names of the two clusters
#' \item ccc_score: communication score calculated on all the interacting ligand/receptor genes in the two clusters
#' \item ngenes_cl1, ngenes_cl2: number of genes in cl1 and cl2, respectively, that participate to the communication score
#' \item nlink: number of links beteewn the genes in cl1 and cl2
#' \item p_value: empirical p-value calculated by using the permutation approcah
#' }
#' }
#' @examples 
#' cl1 <- rep(1, 4)
#' names(cl1) <- c( "A", "B", "C", "D")
#' cl2 <- rep(1, 4)
#' names(cl2) <- c( "E", "F", "G", "H")
#' DEG_list <- list("cl1" = cl1, "cl2" = cl2)
#' adj <- matrix(data = sample(c(0,1), 8*8, replace = T), nrow = 8, ncol = 8, 
#' dimnames = list(LETTERS[1:8], LETTERS[1:8]))
#' ligand <- c("A", "F", "C", "H", "G")
#' receptor <- c("E", "B", "G", "D", "F")
#' cc <- cluster_communication(DEG_list = DEG_list, gene_network_adj = adj, ligand = ligand, 
#' receptor = receptor)
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
  
  perm_list <- parallel::mclapply(1:k, function(x) {
    tmp <- matrix(as.numeric(gene_network_adj), ncol = ncol(gene_network_adj), 
                  dimnames = list(sample(rownames(gene_network_adj), nrow(gene_network_adj))))
    colnames(tmp) <- rownames(tmp)
    return(tmp)
  }, mc.cores = mc_cores_perm
  )
  
  perm_list <- c(list(gene_network_adj), perm_list)
  ans <- parallel::mclapply(1:(k+1), function(x) {
    mat <- perm_list[[x]][names(DEG_list[[1]]), names(DEG_list[[2]]), drop = F]
    matW <- t(DEG_list[[1]]) %*% as.matrix(mat)
    matW <- matW %*% (DEG_list[[2]])
    
    row.col.idx <- which(perm_list[[x]][names(DEG_list[[1]]), names(DEG_list[[2]]), drop = F] == 1, arr.ind = T)
    row.n <- rownames(perm_list[[x]][names(DEG_list[[1]]), names(DEG_list[[2]]), drop = F])[row.col.idx[,1]]
    row.n <- unique(row.n[which(row.n %in% names(DEG_list[[1]]>0))])
    col.n <- colnames(perm_list[[x]][names(DEG_list[[1]]), names(DEG_list[[2]]), drop = F])[row.col.idx[,2]]
    col.n <- unique(col.n[which(col.n %in% names(DEG_list[[2]]>0))])
    
    matW <- data.frame(cl1 = names(DEG_list)[1],
                       cl2 = names(DEG_list)[2],
                       ccc_score = matW,
                       ngenes_cl1 = length(row.n),
                       ngenes_cl2 = length(col.n),
                       nlink = sum(mat),
                       stringsAsFactors = F)
    
    return(matW)
  }, mc.cores = mc_cores_perm)
  
  
  p_list <- parallel::mclapply(ans, function(x) {
    rownames(x) <- paste(x$cl1, x$cl2, sep = "|")
    x <- x[, "ccc_score", drop = F]
    return(x)
  }, mc.cores = mc_cores_perm)
  p_val <- calc_p(p_list)
  
  out <- ans[[1]]
  out$p_value <- p_val
  
  mat <- gene_network_adj[names(DEG_list[[1]]), names(DEG_list[[2]]), drop=F]
  idx <- which(mat==1, arr.ind = T)
  cl_ct <- data.frame(cl1 = rep(names(DEG_list)[1], nrow(idx)),
                      cl1_gene = rownames(idx),
                      cl1_type = rep("ligand", nrow(idx)),
                      cl2 = rep(names(DEG_list)[2], nrow(idx)),
                      cl2_gene = colnames(mat)[idx[,2]],
                      cl2_type = rep("ligand", nrow(idx)),
                      score = DEG_list[[1]][rownames(idx)] * DEG_list[[2]][colnames(mat)[idx[,2]]],
                      stringsAsFactors = F)
  cl_ct$cl1_type[which(cl_ct$cl1_gene %in% receptor)] <- "receptor"
  cl_ct$cl2_type[which(cl_ct$cl2_gene %in% receptor)] <- "receptor"
  cl_ct$cl1_type[which(cl_ct$cl1_gene %in% mixed)] <- "ligand/receptor"
  cl_ct$cl2_type[which(cl_ct$cl2_gene %in% mixed)] <- "ligand/receptor"
  
  
  return(list(communications_info = cl_ct, cc_communications= out))
}
