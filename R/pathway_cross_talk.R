#' Caluclate cross-talks between pathways
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param wight an ordered weight vertex vector
#' @param cluster numebr of threads to be used
#' @import parallel
#' @import doParallel
#' @export

pathway_cross_talk <- function (pathway_list, gene_network_adj, weight = attr_v$Sp, cluster = 4) {
  
  gene_network_adj <- data.frame(sign(gene_network_adj), stringsAsFactors = F)
  colnames(gene_network_adj) <- rownames(gene_network_adj)
  
  comb_p <- expand.grid(ptw1 = names(pathway_list),
                        ptw2 = names(pathway_list))
  comb_p <- data.frame(comb_p[ !duplicated(apply(comb_p, 1, sort), MARGIN = 2), ], 
                       stringsAsFactors = F)
  comb_p <- comb_p[-which(comb_p$ptw1 == comb_p$ptw2),]
  
  
  xx <- apply(comb_p, 1, function (x) {
    idx_1 <-  which(rownames(gene_network_adj) %in% pathway_list[[x["ptw1"]]])
    idx_2 <-  which(rownames(gene_network_adj) %in% pathway_list[[x["ptw2"]]])
    return( gene_network_adj[idx_1, idx_2])
  }
  )
  names(xx) <- paste(comb_p[,1], comb_p[,2], sep = "|")
  idx <- which(lapply(xx, function(n) sum(n)) !=0)
  xx2 <- xx[idx]
  
  pct <- mclapply(names(xx2), function(x) {
    p_i <- unlist(lapply(strsplit(x, "|", fixed = T), "[[", 1))
    p_j <- unlist(lapply(strsplit(x, "|", fixed = T), "[[", 2))
    
    idx_i <- which(rownames(gene_network_adj) %in% pathway_list[[p_i]])
    idx_j <- which(rownames(gene_network_adj) %in% pathway_list[[p_j]])
    i_genes <- rownames(gene_network_adj)[idx_i]
    j_genes <- rownames(gene_network_adj)[idx_j]
    gene_network_adj_ij <- xx2[[x]] 
    if (length(which(i_genes %in% j_genes))>1) {
      
      s <- i_genes[which(i_genes %in% j_genes)]
      sidx_i <- which(i_genes %in% s)
      sidx_j <- which(j_genes %in% s)
      
      gene_network_adj_ij <- gene_network_adj_ij[-sidx_i, -sidx_j]
      i_genes <- i_genes[-sidx_i]
      j_genes <- j_genes[-sidx_j]
      idx_i <- idx_i[-sidx_i]
      idx_j <- idx_j[-sidx_j]
      
    }
    gene_network_adj_ijW <- t(weight[idx_i]) %*% as.matrix(gene_network_adj_ij)
    gene_network_adj_ijW <- gene_network_adj_ijW %*% weight[idx_j]
    gene_network_adj_ijW <- data.frame(pathway_1 = p_i,
                                       pathway_2 = p_j,
                                       pct = gene_network_adj_ijW,
                                       stringsAsFactors = F)
  },   mc.cores = cluster)
  ans <- do.call(rbind, pct)
  
  
  return(ans)
}
