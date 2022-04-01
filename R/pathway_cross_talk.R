#' Function to caluclate cross-talks between pathways
#' @description Calculates cross-talk between all pathway pairs that show a link between their exclusive genes
#' @details The function defines a subset of the `gene_network_adj` for each pathway pair, 
#'  with the exclusive genes of a pathway in the rows and the exclusive genes of the other in the columns. 
#'  Only the pathway pairs that show a link in these subset will be considered. Then, the CT is calcluated, 
#'  as described in the paper. This approach is applied to the original adjacency matrix and to all its permuted versions.
#'  The permutations are obtained by random sampling genes from the adjacency matrix for each pathway.
#'  The p-value and FDR are calculated on the number of links, without considering the weights of the genes.
#'  We also provide a corrected p-value using BH method.
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj adjacency matrix of the whole gene network considered (can be a sparseMatrix)
#' @param weight an vector of weights for each gene the pathway list. If not provided, the function assigns to each gene
#' a weight of 1
#' @param genes target gene to be used for PCT calculation
#' @param mc_cores_pct numebr of threads to be used for pathway cross talk calculation
#' @param mc_cores_perm number of thread to be used in permutations
#' @param k number of permutation 
#' @return The function returns a table with:
#' \itemize{
#'  \item pathway_1, pathway_2: the pathway pair considered
#'  \item pct: the CT value
#'  \item ngenes_pathway_1, ngenes_pathway_2: number of genes involved in `pathway_1` and `pathway_2`, respectively
#'  \item nlink: number of links between the genes in `pathway_1` and `pathway_2`
#'  \item: weight_pathway1, weight_pathway2: cumulative weights of the genes involved in the 
#'  cross-talk between pathway 1 and 2
#'  \item p_value_link: p-value score calcualted by using the permutation approach
#'  \item FDR_link: empirical FDR
#'  \item p_adj_BH: `p-value_link` adjusted using BH method
#'  \item gene_pathway1, gene_pathway2: gene involevd in the CT in `pathway_1` and `pathway_2`, respectively
#' }
#' @examples  
#'  ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#'  adj <- matrix(data = sample(c(0,1), 20*20, replace = TRUE), nrow = 6, 
#'  ncol = 6, dimnames = list(LETTERS[1:20], LETTERS[1:20]))
#'  pct <- pathway_cross_talk(pathway_list = ptw_list, gene_network_adj = adj, weight = NULL, genes = LETTERS[1:6]
#'   mc_cores_pct = 1, mc_cores_perm = 1, k = 9)
#' @import parallel
#' @import igraph
#' @importFrom stringi stri_c
#' @importFrom reshape2 acast
#' @import Matrix
#' @export
pathway_cross_talk <- function (pathway_list, gene_network_adj, genes, 
                                weight =NULL, 
                                mc_cores_pct = 2, mc_cores_perm = 1, 
                                k = 9) {
  if(is.null(weight) ) {
    weight <- rep(1, length(genes))
    names(weight) <- genes
  } 
  if(!is(gene_network_adj, "sparseMatrix" )) {
    gene_network_adj <- as(gene_network_adj, "dgCMatrix")
  }
  gene_network_adj <- sign(gene_network_adj)
  sub_adj_mt <- gene_network_adj[as.character(genes), as.character(genes), drop = F]
  
  comb_p <- expand.grid(names(pathway_list), names(pathway_list))
  comb_p <- graph_from_edgelist(as.matrix(comb_p), directed = F)
  comb_p <- igraph::simplify(comb_p, remove.multiple = T, remove.loops = T)
  comb_p <- as_edgelist(comb_p)
  
  xx <- mclapply(1:nrow(comb_p),  function (x) {
    
    g.1 <- pathway_list[[comb_p[x,1]]]
    g.2 <- pathway_list[[comb_p[x,2]]]
    idx_1 <- g.1[!g.1 %in% g.2]
    idx_2 <- g.2[!g.2 %in% g.1]
    
    tmp <- sub_adj_mt[idx_1, idx_2, drop = F]
    out <- list(tmp, idx_1, idx_2)
    return( out) 
  }, mc.cores = mc_cores_pct
  )
  
  idx <- which(mclapply(xx, function(n) sum(n[[1]], na.rm = T), mc.cores = mc_cores_pct) !=0)
  
  if(length(idx) == 0) {
    print("no available PCT")
    return("no available PCT")
    
  } else {
    xx <- xx[idx]
    comb_p <- comb_p[idx,, drop = F]
    len1 <- array(unlist(mclapply(xx, function(n) length(n[[2]]), mc.cores = mc_cores_pct)), 
                  dim = c(length(xx), 1))
    len2 <- array(unlist(mclapply(xx, function(n) length(n[[3]]), mc.cores = mc_cores_pct)), 
                  dim = c(length(xx), 1))
    comb_p_len <- matrix(cbind(len1, len2), ncol = 2)
    unq_comb_p_len <- unique(comb_p_len)
    perm_list <- mclapply(1:nrow(unq_comb_p_len), function(x) {
      n.1 <- unq_comb_p_len[x, 1]
      n.2 <- unq_comb_p_len[x, 2]
      out <- perm_link(r = n.1, c = n.2, gene_network_adj, core = mc_cores_perm, k, hash = T)
      return(out)
    }, mc.cores = mc_cores_pct)
    
    unq_comb_p_len <- data.frame(unq_comb_p_len)
    unq_comb_p_len$var <- 1:nrow(unq_comb_p_len)
    unq_comb_p_len <- acast(unq_comb_p_len, X1 ~ X2, value.var = "var")
    
    pct <- mclapply(1:length(xx), function(x) {
      tab <- xx[[x]][[1]]
      n.1 <- length(xx[[x]][[2]])
      n.2 <- length(xx[[x]][[3]])
      nlink <- matrix(sum(tab))
      idx <- unq_comb_p_len[as.character(n.1), as.character(n.2)]
      perm_l <- c(list(nlink), perm_list[[idx]])
      p_val <- calc_p(perm_l)
      ct <- cross_talk(mat = tab, weight = weight)
      out <- array(c(comb_p[x, 1], comb_p[x, 2], ct[1:6], as.vector(p_val), ct[7:8]), dim = c(1, 11))
      return(out)
    },   mc.cores = mc_cores_pct)
    
    pct <- do.call(rbind, pct)
    all.v <- unlist(perm_list)
    all.v <- c(as.numeric(pct[, 6]), all.v)
    link_FDR <- eFDR(real_values = as.numeric(pct[, 6]), all_values = all.v, mc.cores = mc_cores_pct)
    p_bh <- stats::p.adjust(pct[, 9], method = "BH")
    pct <- cbind(pct[,1:9], link_FDR, p_bh,  pct[, 10:11])
    colnames(pct) <- c("pathway_1", "pathway_2", "pct", "ngenes_pathway1",
                       "ngenes_pathway2", "nlink","weight_pathway1", 
                       "weight_pathway2", "p_value_link", "FDR_link", "p_adj_BH", "gene_pathway1", "gene_pathway2")
    pct <- data.frame(pct, stringsAsFactors = F)
    return(pct)
    
    
  }
  
  
}
