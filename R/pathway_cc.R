#' Caluclate cross-talks between pathways
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param wight an ordered weight vertex vector
#' @param cluster numebr of threads to be used
#' @import parallel
#' @import doParallel
#' @import igraph
#' @export


pathway_cc <- function (pathway_list, gene_network_adj, weight = attr_v$Sp, cluster = 4) {
  
  gene_network_adj <- data.frame(sign(gene_network_adj), stringsAsFactors = F)
  colnames(gene_network_adj) <- rownames(gene_network_adj)
  
  xx <- lapply(pathway_list, function(x) {
    idx <- which(rownames(gene_network_adj) %in% x)
    return(gene_network_adj[idx,idx])
  })
  
  idx <- which(lapply(xx, function(n) sum(n)) !=0)
  xx2 <- xx[idx]
  
  pct <- mclapply(names(xx2), function(x) {
    idx <- which(rownames(gene_network_adj) %in% pathway_list[[x]])
    names(idx) <- rownames(gene_network_adj)[idx]
    gene_network_adj_ij <- gene_network_adj[idx, idx]
    cc <- igraph::components(graph = igraph::graph_from_adjacency_matrix(as.matrix(gene_network_adj_ij),
                                                                         mode = "undirected"))
    cc_list <- list(components_results = cc,
                    pathway_cc = matrix(data=0, nrow = 1, 
                                        ncol = cc$no, dimnames = list("value",
                                                                      paste("CC", 1:cc$no, sep = "_"))))
    
    for(j in 1:cc$no) {
      cc_genes <- names(which(cc$membership == j))
      CCgene_network_adj_ij <- gene_network_adj_ij[cc_genes, cc_genes]
      gene_network_adj_ijW <- sum(t(weight[idx[cc_genes]]) %*% as.matrix(CCgene_network_adj_ij))
      gene_network_adj_ijW <- gene_network_adj_ijW/length(cc_genes)
      cc_list$pathway_cc[1,j] <- gene_network_adj_ijW
      
    }
    return(cc_list)
    
  },   mc.cores = cluster)
  
  names(pct) <- names(xx2)
  
  return(pct)
}
