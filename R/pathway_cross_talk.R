#' Pathway Cross Talk
#' @param pathway_list list of gene sets
#' @param gene_network_adj adjacency matrix of the gene network under analysis
#' @export

pathway_cross_talk <- function (pathway_list, gene_network_adj) {
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
      if (i < j) {
        if (length(which(i_genes %in% j_genes))>1) {
          s <- i_genes[which(i_genes %in% j_genes)]
          idx_i <- which(i_genes %in% s)
          idx_j <- which(j_genes %in% s)
          #sh <- rowSums(matrix(gene_network_adj_ij[idx_i, -idx_j])) + colSums(matrix(gene_network_adj_ij[-idx_i, idx_j]))
          
          gene_network_adj_ij <- gene_network_adj_ij[-idx_i, -idx_j]
          
          #tmp <- sum(gene_network_adj_ij) + length(which(sh != 0))
          
        }
        #print(i)
        #print(j)
        n <- sum(gene_network_adj_ij)
        #nab <- factorial((length(i_genes) + length(j_genes)))/(factorial(2)* factorial(abs(length(i_genes) + length(j_genes))-2))
        #na <- factorial(length(i_genes))/(factorial(2)* factorial(abs(length(i_genes)-2)))
        #nb <- factorial(length(j_genes))/(factorial(2)* factorial(abs(length(j_genes)-2)))
        d <- choose(length(i_genes) + length(j_genes), 2) - (choose(length(i_genes), 2) + choose(length(j_genes),2))
          #nab - (na + nb)
        ans[i, j] <- n/d
      }
    }
  }
  return(ans)
}

