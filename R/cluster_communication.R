#' Function to caluclate communication between cell clusters
#' @description cluster_communication calculates a score for the communication between all the possible pairs 
#' in the cluster list provided
#' @details The CT formula is used to calculate a score for the communication between all possible cluster pairs. 
#' The same approach is applied also to `k` permuted version of each cluster, obtained by random sampling of genes
#' in `gene_network_adj`. The p-value and FDR are the calculated by comparing the number of links between each 
#' cluster pairs and their permuted versions
#' @param cl_list the output of `preparing_cl_list()` or a gene list composed by a vector for each cluster where each vector is
#' a named gene rank or weight vector.
#' @param gene_network_adj adjacency matrix of the whole gene network
#' @param mc_cores_perm number of threads to be used to calcualte permutation
#' @param mc_cores_ccc number of threads to be used for CCC calculation
#' @param k number of permutation
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
#' \item p_value_link: empirical p-value calculated by using the permutation approcah
#' \item FDR_link: empirical false discovery rate
#' \item p_adj_BH: `p-value_link` correctd by using BH method
#' \item gene_weight_cl1, gene_weight_cl2: cumulative weights of the genes involved in the 
#' communication in cluster 1 and 2
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
#' @importFrom stringi stri_c
#' @import Matrix
#' @export

cluster_communication <- function(cl_list, gene_network_adj, k = 9, 
                                  mc_cores_perm = 1, mc_cores_ccc = 1) {
  
  if(!is(gene_network_adj, "sparseMatrix" )) {
    gene_network_adj <- as(gene_network_adj, "dgCMatrix")
  }
  gene_network_adj <- sign(gene_network_adj)
  
  comb_p <- expand.grid(names(cl_list), names(cl_list))
  g <- graph_from_edgelist(as.matrix(comb_p), directed = F)
  g <- igraph::simplify(g, remove.multiple = T, remove.loops = T)
  comb_p <- as_edgelist(g)
  l_cl_list <- lengths(cl_list)
  names(l_cl_list) <- names(cl_list)
  
  perm_list <- mclapply(1:nrow(comb_p), function(x) {
    n.1 <- l_cl_list[comb_p[x,1]]
    n.2 <- l_cl_list[comb_p[x,1]]
    out <- perm_link(r = n.1, c = n.2, gene_network_adj, core = mc_cores_perm, k, hash = F)
    return(out)
  }, mc.cores = mc_cores_ccc)
  
  
  ans <- parallel::mclapply(1:nrow(comb_p), function(z) {
    g.1 <- cl_list[[comb_p[z, 1]]]
    n.g1 <- as.character(names(g.1))
    n.1 <- length(g.1)
    g.2 <- cl_list[[comb_p[z, 2]]]
    n.g2 <- as.character(names(g.2))
    n.2 <- length(g.2)
    
    tab <- gene_network_adj[n.g1, n.g2]
    nlink <- matrix(sum(tab))
    perm_l <- perm_list[[z]]
    all_perm <- c(list(nlink), perm_l)
    p_val <- calc_p(all_perm)
    ct <- cross_talk(mat = tab, weight = c(g.1, g.2))
    out <- array(c(comb_p[z, 1], comb_p[z, 2], ct[1:4], as.vector(p_val), ct[5:8]), dim = c(1, 11))
    
    return(out)
  }, mc.cores = mc_cores_perm)
  
  ans <- do.call(rbind, ans)
  all.v <- unlist(perm_list)
  all.v <- c(ans[,6], all.v)
  link_FDR <- eFDR(real_values = ans[,6], all_values = all.v, mc.cores = mc_cores_ccc)
  p_val_BH <- stats::p.adjust(ans[,7], method = "BH")
  ans <- cbind(ans[, 1:7], link_FDR, p_val_BH, ans[,8:11])
  
  colnames(ans) <- c("cl1", "cl2", "ccc_score", "ngenes_cl1", "ngenes_cl2",
                     "nlink","p_value_link", "link_FDR", "p_adj_BH", "gene_weight_cl1", 
                     "gene_weight_cl2", "genes_cl1", "genes_cl2")
  ans <- data.frame(cl1 = ans[, 1], 
                    cl2 = ans[, 2], 
                    ccc_score = as.numeric(ans[, 3]), 
                    ngenes_cl1 = as.numeric(ans[, 4]), 
                    ngenes_cl2 = as.numeric(ans[, 5]),
                    nlink = as.numeric(ans[, 6]), 
                    p_value_link = as.numeric(ans[, 7]), 
                    link_FDR = as.numeric(ans[, 8]), 
                    p_adj_BH = as.numeric(ans[, 9]), 
                    gene_weight_cl1 = as.numeric(ans[, 10]), 
                    gene_weight_cl2 = as.numeric(ans[, 11]), 
                    genes_cl1 = ans[, 12], 
                    genes_cl2 = ans[, 13], 
                    stringsAsFactors = F)
  
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
  ct_info <- data.frame(cl1 = ct_info[, 1],
                        cl1_gene = ct_info[, 2],
                        cl2 = ct_info[, 3],
                        cl2_gene = ct_info[, 4],
                        score = as.numeric(ct_info[, 5]),
                        stringsAsFactors = F)
  
  return(list(communications_info = ct_info, cc_communications= ans))
}



