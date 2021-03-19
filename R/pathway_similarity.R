#' Pathway similarity
#' @param pathway_list list of gene sets
#' @param gene_network_adj adjacency matrix of the gene network under analysis
#' @param sim_method method to define the similarity between pathways
#' @export

pathway_similarity <- function(pathway_list, gene_network_adj, sim_method = "overlap") {
  ans <- matrix(0, nrow = length(pathway_list), ncol = length(pathway_list), 
                dimnames = list(names(pathway_list), names(pathway_list)))
  gene_network_adj <- sign(gene_network_adj)
  for (i in 1:nrow(ans)) {
    i_genes <- rownames(gene_network_adj)[rownames(gene_network_adj) %in% 
                                            pathway_list[[i]]]
    for (j in 1:ncol(ans)) {
      j_genes <- rownames(gene_network_adj)[rownames(gene_network_adj) %in% 
                                              pathway_list[[j]]]
      gene_network_adj_ij <- gene_network_adj[rownames(gene_network_adj) %in% 
                                                i_genes, rownames(gene_network_adj) %in% j_genes]
      if (i > j) {
        if (any(i_genes %in% j_genes)) {
          ans[i, j] <- calc_set_similarity(pathway_list[[i]], pathway_list[[j]], method = sim_method)
        }
      }
    }
  }
  return(ans)
}