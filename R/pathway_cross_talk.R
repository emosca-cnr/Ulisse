#' Caluclate cross-talks between pathways
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param wight an ordered weight vertex vector
#' @param mc_cores_pct numebr of threads to be used for pathway cross talk calculation
#' @param mc_cores_permutation number of thread to be used in permutations
#' @param k number of permutation of the adjacency matrix
#' @import parallel
#' @import igraph
#' @export
pathway_cross_talk <- function (pathway_list, gene_network_adj, 
                                weight = attr_v$Sp, 
                                mc_cores_pct = 2, mc_cores_perm = 1, #cluster = 4, 
                                k = 9) {
  gene_network_adj <- sign(gene_network_adj)
  names(weight) <- rownames(gene_network_adj)
  
  comb_p <- expand.grid(names(pathway_list), names(pathway_list))
  g <- graph_from_edgelist(as.matrix(comb_p), directed = F)
  g <- igraph::simplify(g, remove.multiple = T, remove.loops = T)
  comb_p <- as_edgelist(g)
  xx <- mclapply(1:nrow(comb_p),  function (x) {
    
    idx_1 <- as.character(pathway_list[[comb_p[x,1]]][!pathway_list[[comb_p[x,1]]] %in% pathway_list[[comb_p[x,2]]]])
    idx_2 <- as.character(pathway_list[[comb_p[x,2]]][!pathway_list[[comb_p[x,2]]] %in% pathway_list[[comb_p[x,1]]]])
    
    tmp <- gene_network_adj[idx_1, idx_2, drop = F]
    
    return( tmp) 
  }, mc.cores = mc_cores_pct
  )
  
  idx <- which(mclapply(xx, function(n) sum(n, na.rm = T), mc.cores = mc_cores_pct) !=0)
  
  if(length(idx) == 0) {
    print("no available PCT")
    return("no available PCT")
  } else {
    xx2 <- xx[idx]
    comb_p <- comb_p[idx,]
    
    perm_list <- mclapply(1:k, function(x) {
      tmp <- matrix(as.numeric(gene_network_adj), ncol = ncol(gene_network_adj), 
                    dimnames = list(sample(rownames(gene_network_adj), nrow(gene_network_adj))))
      colnames(tmp) <- rownames(tmp)
      return(tmp)
    }, mc.cores = mc_cores_perm
    )
    xxPL <- mclapply(1:k, function(j) mclapply(1:nrow(comb_p),  function (x) {
      
      idx_1 <- as.character(pathway_list[[comb_p[x,1]]][!pathway_list[[comb_p[x,1]]] %in% pathway_list[[comb_p[x,2]]]])
      idx_2 <- as.character(pathway_list[[comb_p[x,2]]][!pathway_list[[comb_p[x,2]]] %in% pathway_list[[comb_p[x,1]]]])
      
      tmp <- perm_list[[j]][idx_1, idx_2, drop = F]
      return( tmp) 
    }, mc.cores = mc_cores_pct
    ), mc.cores = mc_cores_perm)
    all_pct <- c(list(xx2), xxPL)
    
    pct <- mclapply(1:(k+1), function(j) {
      tmp <-  mclapply(1:length(all_pct[[j]]), function(x) {
        
        gene_network_adj_ijW <- t(weight[rownames(all_pct[[j]][[x]])]) %*% as.matrix(all_pct[[j]][[x]])
        gene_network_adj_ijW <- gene_network_adj_ijW %*% weight[colnames(all_pct[[j]][[x]])]
        gene_network_adj_ijW <- data.frame(pathway_1 = comb_p[x,1],
                                           pathway_2 = comb_p[x,2],
                                           pct = gene_network_adj_ijW,
                                           ngenes_pathway1 = length(pathway_list[[comb_p[x,1]]]),
                                           ngenes_pathway2 = length(pathway_list[[comb_p[x,2]]]),
                                           nlink = sum(all_pct[[j]][[x]]),
                                           stringsAsFactors = F)
      },   mc.cores = mc_cores_pct)
      ans <- do.call(rbind, tmp)
      return(ans)
    }, mc.cores = mc_cores_perm
    )
    p_list <- mclapply(pct, function(x) {
      rownames(x) <- paste(x$pathway_1, x$pathway_2, sep = "|")
      x <- x[, "pct", drop = F]
      return(x)
    }, mc.cores = mc_cores_perm)
    
    p_val <- calc_p(p_list)
    
    out <- pct[[1]]
    out$p_value <- p_val
    out$gene_pathway1 <- rep(NA, nrow(out))
    out$gene_pathway2 <- rep(NA, nrow(out))
    for (i in 1:nrow(out)) {
      out$gene_pathway1[i] <- paste(unlist(pathway_list[out$pathway_1[i]]), collapse = ";")
      out$gene_pathway2[i] <- paste(unlist(pathway_list[out$pathway_2[i]]), collapse = ";")
    }
    out$eFDR <- eFDR(real_values = as.vector(unlist(p_list[[1]])), all_values = as.vector(unlist(p_list)))
    return(out)
    
    
  }
  
  
}
