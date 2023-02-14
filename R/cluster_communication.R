#' Function to calculate communication between cell clusters
#' @description `cluster_communication()` calculates a communication score between all the possible pairs 
#' in the cluster list provided
#' @details The function calculate a communication score between all possible cluster pairs that can be obtained from the `cl_list` provided by using the 
#' weights or rank provided in the lists. For each cluster pair, the function samples randomly from the `gene_network_adj` two lists of the same dimensions of the 
#' lists of the cluster pair considered, and calculates the number of links between them. This process is repeated `k` times. The number of links obtained with the 
#' `k` permuted version of the cluster pair lists are compared to the real number of links between the clusters to calculate the p-value and FDR. 
#' @param cl_list the output of `preparing_cl_list()` or a gene list composed by a named vector. Each vector should be composed by the ranks or weights 
#' named after the respective genes. The element of the list should be named after the name of the clusters. These name will be then used in the output.
#' If no weights are available use 1 for each gene. In this case the communication score will be equal the number of links 
#' between the clusters
#' @param gene_network_adj adjacency matrix of the whole communication network
#' @param mc_cores_perm number of threads to be used to calculate permutations
#' @param mc_cores_ccc number of threads to be used for CCC calculation
#' @param k number of permutations
#' @return The function returns a list of two objects:
#' \enumerate{
#' \item communications_info: a data.frame with the detailed description of the communication between each cluster pairs
#' \itemize{
#' \item cl1, cl2: name of the clusters (derived from the name of the list in `cl_list`)
#' \item cl1_gene, cl2_gene: genes involved in the communication
#' \item score: the score of the communication between cl1_gene and cl2_gene, calculated considering their weights (if provided)
#' }
#' \item {cc_communications}: a data.frame with the cumulative communication score
#' \itemize{
#' \item cl1, cl2: names of the two clusters
#' \item ccc_score: communication score calculated on all the interacting genes in the two clusters. This value is equal to the sum of all 
#' the `communications_info$score` for a cluster pair
#' \item ngenes_cl1, ngenes_cl2: number of genes in `cl1` and `cl2`, respectively, involved in the communication score
#' \item nlink: number of links between the genes in `cl1` and `cl2`
#' \item p_value_link: empirical p-value calculated by using the permutation approach
#' \item FDR_link: empirical false discovery rate
#' \item p_adj_BH: `p-value_link` corrected by using BH method (see `p.adjust()` function)
#' \item gene_weight_cl1, gene_weight_cl2: cumulative weights of the genes involved in the 
#' communication in `cl1` and `cl2`
#' \item genes_cl1, genes_cl2: list of the genes involved in the communication between `cl1` and `cl2`
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
  
  sub_adj <- mclapply(1:nrow(comb_p), function(x) {
    g.1 <- cl_list[[comb_p[x, 1]]]
    n.g1 <- as.character(names(g.1))
    n.1 <- length(g.1)
    g.2 <- cl_list[[comb_p[x, 2]]]
    n.g2 <- as.character(names(g.2))
    n.2 <- length(g.2)
    
    tab <- gene_network_adj[n.g1, n.g2]
    return(tab)
  })
  
  idx <- which(mclapply(sub_adj, function(n) sum(n, na.rm = T), mc.cores = mc_cores_ccc) !=0)
  
  if(length(idx) == 0 ) {
    print("No CCC available")
    return("No CCC available")
  } else {
    
    sub_adj <- sub_adj[idx]
    comb_p <- comb_p[idx,]
    perm_list <- mclapply(1:nrow(comb_p), function(x) {
      n.1 <- l_cl_list[comb_p[x,1]]
      n.2 <- l_cl_list[comb_p[x,2]]
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
      ct <- cross_talk(mat = tab, weight = list(g1 = g.1, g2 =g.2))
      out <- array(c(comb_p[z, 1], comb_p[z, 2], ct[1:4], as.vector(p_val), ct[5:8]), dim = c(1, 11))
      
      return(out)
    }, mc.cores = mc_cores_perm)
    
    ans <- do.call(rbind, ans)
    all.v <- unlist(perm_list)
    all.v <- c(as.numeric(ans[,6]), all.v)
    link_FDR <- eFDR(real_values = as.numeric(ans[,6]), all_values = all.v, mc.cores = mc_cores_ccc)
    p_val_BH <- stats::p.adjust(as.numeric(ans[,7]), method = "BH")
    ans <- data.frame(cl1 = ans[, 1], 
                      cl2 = ans[, 2], 
                      ccc_score = as.numeric(ans[, 3]), 
                      ngenes_cl1 = as.numeric(ans[, 4]), 
                      ngenes_cl2 = as.numeric(ans[, 5]),
                      nlink = as.numeric(ans[, 6]), 
                      p_value_link = as.numeric(ans[, 7]), 
                      FDR_link = link_FDR, 
                      p_adj_BH = p_val_BH, 
                      gene_weight_cl1 = as.numeric(ans[, 8]), 
                      gene_weight_cl2 = as.numeric(ans[, 9]), 
                      genes_cl1 = ans[, 10], 
                      genes_cl2 = ans[, 11], 
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
  
  
}



