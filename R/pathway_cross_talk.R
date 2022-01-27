#' Function to caluclate cross-talks between pathways
#' @description Calculates cross-talk between all pathway pairs that show a link between their exclusive genes
#' @details The function defines a subset of the `gene_network_adj` for each pathway pair, 
#'  with the exclusive genes of a pathway in the rows and the exclusive genes of the other in the columns. 
#'  Only the pathway pairs that show a link in these subset will be considered. Then, the CT is calcluated, 
#'  as described in the paper. This approach is applied to the original adjacency matrix and to all its permuted versions.
#'  The CT values and the permuted CT are then used to calculated the empirical p-value and the FDR score
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param weight an vector of weights for each gene in gene_network_adj
#' @param mc_cores_pct numebr of threads to be used for pathway cross talk calculation
#' @param mc_cores_perm number of thread to be used in permutations
#' @param k number of permutation of the adjacency matrix
#' @return The function returns a table with:
#' \itemize{
#'  \item pathway_1, pathway_2: the pathway pair considered
#'  \item pct: the CT value
#'  \item ngenes_pathway_1, ngenes_pathway_2: number of genes involved in `pathway_1` and `pathway_2`, respectively
#'  \item nlink: number of links between the genes in `pathway_1` and `pathway_2`
#'  \item gene_pathway1, gene_pathway2: gene involevd in the CT in `pathway_1` and `pathway_2`, respectively
#'  \item p_value: p-value score calcualted by using the permutation approach
#'  \item eFDR: empirical FDR
#' }
#' @examples  
#'  ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#'  adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#'  ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#'  wgt <- rep(1, 6)
#'  pct <- pathway_cross_talk(pathway_list = ptw_list, gene_network_adj = adj, weight = wgt, 
#'   mc_cores_pct = 1, mc_cores_perm = 1, k = 9)
#' @import parallel
#' @import igraph
#' @export
pathway_cross_talk <- function (pathway_list, gene_network_adj, 
                                weight, 
                                mc_cores_pct = 2, mc_cores_perm = 1, 
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
        row.col.idx <- which(all_pct[[j]][[x]] == 1, arr.ind = T)
        row.n <- rownames(all_pct[[j]][[x]])[row.col.idx[,1]]
        row.n <- unique(row.n[which(row.n %in% names(weight>0))])
        col.n <- colnames(all_pct[[j]][[x]])[row.col.idx[,2]]
        col.n <- unique(col.n[which(col.n %in% names(weight>0))])
        
        gene_network_adj_ijW <- data.frame(pathway_1 = comb_p[x,1],
                                           pathway_2 = comb_p[x,2],
                                           pct = gene_network_adj_ijW,
                                           ngenes_pathway1 = length(row.n),
                                           ngenes_pathway2 = length(col.n),
                                           nlink = sum(all_pct[[j]][[x]][which(weight[rownames(all_pct[[j]][[x]])] >0), which(weight[colnames(all_pct[[j]][[x]])] > 0)]),
                                           gene_pathway1 = paste(row.n, collapse = ";"),
                                           gene_pathway2 = paste(col.n, collapse = ";"),
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
    p_val <- as.vector(calc_p(p_list))
    
    out <- pct[[1]]
    out$p_value <- p_val
    
    out$eFDR <- eFDR(real_values = as.vector(unlist(p_list[[1]])), all_values = as.vector(unlist(p_list)))
    return(out)
    
    
  }
  
  
}
