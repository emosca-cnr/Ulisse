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
#' @import stringi
#' @export
pathway_cross_talk <- function (pathway_list, gene_network_adj, 
                                weight, 
                                mc_cores_pct = 2, mc_cores_perm = 1, 
                                k = 9) {
  gene_network_adj <- sign(gene_network_adj)
  names(weight) <- rownames(gene_network_adj)
  weight <- weight[weight !=0]
  gene_network_adj[names(weight), names(weight)]
  
  comb_p <- expand.grid(names(pathway_list), names(pathway_list))
  comb_p <- graph_from_edgelist(as.matrix(comb_p), directed = F)
  comb_p <- igraph::simplify(comb_p, remove.multiple = T, remove.loops = T)
  comb_p <- as_edgelist(comb_p)
  
  xx <- mclapply(1:nrow(comb_p),  function (x) {
    
    g.1 <- pathway_list[[comb_p[x,1]]]
    g.2 <- pathway_list[[comb_p[x,2]]]
    idx_1 <- as.character(g.1[!g.1 %in% g.2])
    idx_2 <- as.character(g.2[!g.2 %in% g.1])
    
    tmp <- gene_network_adj[idx_1, idx_2, drop = F]
    return( tmp) 
  }, mc.cores = mc_cores_pct
  )
  
  idx <- which(mclapply(xx, function(n) sum(n, na.rm = T), mc.cores = mc_cores_pct) !=0)
  
  if(length(idx) == 0) {
    print("no available PCT")
    return("no available PCT")
    
  } else {
    xx <- xx[idx]
    comb_p <- comb_p[idx,]
    
    perm_list <- mclapply(1:k, function(x) {
      tmp <- gene_network_adj
      rownames(tmp) <- sample(rownames(tmp), nrow(tmp))
      colnames(tmp) <- rownames(tmp)
      
      xxPL <- mclapply(1:length(xx),  function (j) {
        g.1 <- rownames(xx[[j]])
        g.2 <- colnames(xx[[j]])
        
        tmp <- tmp[g.1, g.2, drop = F]
        
        return( tmp) 
      }, mc.cores = mc_cores_pct)
      return(xxPL)
    }, mc.cores = mc_cores_perm
    )
    
    all_pct <- c(list(xx), perm_list)
    rm(xx, perm_list)
    
    pct <- mclapply(1:(k+1), function(j) {
      tmp <-  mclapply(1:length(all_pct[[j]]), function(x) {
        g.1 <- rownames(tmp)
        g.2 <- colnames(tmp)
        gene_network_adj_ijW <- t(weight[g.1]) %*% as.matrix(tmp)
        gene_network_adj_ijW <- gene_network_adj_ijW %*% weight[g.2, drop = F]
        row.col.idx <- which(tmp == 1, arr.ind = T)
        row.n <- g.1[row.col.idx[,1]]
        row.n <- row.n[row.n %in% names(weight)]
        col.n <- g.2[row.col.idx[,2]]
        col.n <- col.n[col.n %in% names(weight)]
        
        #mat.out <- matrix(data = NA, nrow = 1, ncol = 8)
        mat.out <- array(data = NA, dim = 8)
        mat.out[1] = comb_p[x,1]
        mat.out[2] = comb_p[x,2]
        mat.out[3] = gene_network_adj_ijW
        mat.out[4] = length(row.n)
        mat.out[5] = length(col.n)
        mat.out[6] = sum(tmp)
        mat.out[7] = stri_c(row.n, collapse = ";")
        mat.out[8] = stri_c(col.n, collapse = ";")
        
        return(mat.out)
      },   mc.cores = mc_cores_pct)
      ans <- do.call(rbind, tmp)
      colnames(ans) <- c("pathway_1", "pathway_2", "pct", "ngenes_pathway1",
                         "ngenes_pathway2", "nlink","weight_pathway1", 
                         "weight_pathway2", "gene_pathway1", "gene_pathway2")
      return(ans)
    }, mc.cores = mc_cores_perm
    )
    
    p_list <- mclapply(pct, function(x) {
      tmp <- matrix(as.numeric(x[,"pct"]), ncol = 1)
      rownames(tmp) <- paste(x[,1], x[,2], sep = "|")
      return(tmp)
    }, mc.cores = mc_cores_perm)
    p_val <- calc_p(p_list)
    
    out <- data.frame(pct[[1]], stringsAsFactors = F)
    out$p_value <- as.vector(p_val)
    out$eFDR <- eFDR(real_values = as.vector(unlist(p_list[[1]])), all_values = as.vector(unlist(p_list)))
    return(out)
    
    
  }
  
  
}
