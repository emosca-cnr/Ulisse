#' Cross-talks between gene-sets
#' @description The function calculates cross-talk between all gene-set pairs that show a link between their genes
#' @details The function takes as inputs the adjacency matrix of the biological network (`gene_network_adj`) and the gene-set list, 
#'  either obtained through `preparing_gs_list`, `preparing_msigdb_list()`, `preparing_DEG_list()` or `preparing_expr_list()` functions. 
#'  For each gene-set pair that shows at least a link between them, the function uses the adjacency matrix to calculate cross-talk as described in the paper. 
#'  Then, for each  pair the function samples randomly from the `gene_network_adj` two lists of the same dimensions of the gene-sets and calculates the 
#'  number of links between them. This process is repeated `k` times. The number of links obtained with the `k` permuted version of the gene-set pairs 
#'  are compared to the real to calculate the p-value and FDR. We also provide a corrected p-value using BH method.
#' @param gs_list a named list of genes grouped into gene sets, as obtained from `preparing_gs_list`, `preparing_msigdb_list()`, 
#'  `preparing_DEG_list()` or `preparing_expr_list()` functions
#' @param gene_network_adj adjacency matrix of the whole gene network considered (can be a sparseMatrix)
#' @param k number of permutations
#' @param shared logical, if the cross-talk calculation should consider shared genes (`TRUE`, suggested only for cell-cell communication) or 
#'  not (`FALSE`, suggested for pathway cross-talk)
#' @param hash logical, if hash map should be used to speed cross-talk calculation. Suggested only with high number of gene-sets (like in pathway cross-talk)
#' @param ct_info logical, if the function has to return the detail of gene-gene interactions between gene sets
#' @param mc_cores_ct number of threads to be used for cross talk calculation
#' @param mc_cores_perm number of thread to be used in permutations
#' @return The function can return two output type:
#' \itemize{
#'  \item `ct_info = FALSE`: only the cross-talk results are returned as output
#'  \itemize{
#'    \item gs1,gs2: names of gene-set pairs considered
#'    \item ct_score: cross-talk score
#'    \item ngenes_gs1,ngenes_gs2: number of genes in `gs1` and `gs2`
#'    \item nlink: number of links between `gs1` and `gs2`
#'    \item p_value_link: p-value score calculated by using the permutation approach
#'    \item FDR_link: empirical FDR
#'    \item p_adj_BH: `p-value_link` adjusted using BH method (see `p.adjust()` function)
#'    \item weight_gs1,weight_gs2: cumulative weights of the genes involved in cross-talk between `gs1` and `gs2`
#'    \item genes_gs1,genes_gs2: genes involved in cross-talk between `gs1` and `gs2`
#'  }
#'  \item `ct_info = TRUE`: the cross-talk results and the detail of gene-gene interactions between gene-sets are returned as output
#'  \enumerate{
#'   \item ct_info: gene-gene interaction details
#'   \itemize{
#'    \item gs1,gs2: name of the gene-sets
#'    \item gs1_gene,gs2_gene: genes involved in the cross-talk
#'    \item score: the score of the interaction between `gs1_gene` and `gs1_gene`, calculated by multiplying their weigths
#'   }
#'  \item ct_res: cross-talk results
#'  \itemize{
#'    \item gs1,gs2: names of gene-set pairs considered
#'    \item ct_score: cross-talk score
#'    \item ngenes_gs1,ngenes_gs2: number of genes in `gs1` and `gs2`
#'    \item nlink: number of links between `gs1` and `gs2`
#'    \item p_value_link: p-value score calculated by using the permutation approach
#'    \item FDR_link: empirical FDR
#'    \item p_adj_BH: `p-value_link` adjusted using BH method (see `p.adjust()` function)
#'    \item weight_gs1,weight_gs2: cumulative weights of the genes involved in cross-talk between `gs1` and `gs2`
#'    \item genes_gs1,genes_gs2: genes involved in cross-talk between `gs1` and `gs2`
#'    }
#'  }
#' }
#' @examples
#' \dontrun{  
#'  ptw_list <- list(ptwA = setNames(c(1, 1, 1), c("A", "B","C")), 
#'  ptwB = setNames(c(1, 1, 1), c("D", "E", "F")), ptwC = setNames(c(1, 1, 1), c("A", "B", "E")))
#'  adj <- matrix(data = sample(c(0,1), 20*20, replace = TRUE), nrow = 6, 
#'  ncol = 6, dimnames = list(LETTERS[1:20], LETTERS[1:20]))
#'  pct <- gs_cross_talk(pathway_list = ptw_list, gene_network_adj = adj, 
#'                       shared = FALSE, hash = FALSE, ct_info = FALSE, 
#'                       mc_cores_pct = 1, mc_cores_perm = 1, k = 9)
#'  }
#' @import parallel
#' @import igraph
#' @importFrom stringi stri_c
#' @importFrom reshape2 acast
#' @import Matrix
#' @export


gs_cross_talk <- function(gs_list, gene_network_adj, k = 9, shared = F, hash = T, ct_info = F,
                                  mc_cores_perm = 1, mc_cores_ct = 1) {
  
  if(!is(gene_network_adj, "sparseMatrix" )) {
    gene_network_adj <- as(gene_network_adj, "dgCMatrix")
  }
  gene_network_adj <- sign(gene_network_adj)
  gs_list <- lapply(gs_list, function(x) x <- x[names(x) %in% rownames(gene_network_adj)])
  
  comb_p <- expand.grid(names(gs_list), names(gs_list))
  g <- graph_from_edgelist(as.matrix(comb_p), directed = F)
  g <- igraph::simplify(g, remove.multiple = T, remove.loops = T)
  comb_p <- as_edgelist(g)
  l_cl_list <- lengths(gs_list)
  names(l_cl_list) <- names(gs_list)
  
  sub_adj <- mclapply(1:nrow(comb_p), function(x) {
    g.1 <- gs_list[[comb_p[x, 1]]]
    n.g1 <- as.character(names(g.1))
    #n.1 <- length(g.1)
    g.2 <- gs_list[[comb_p[x, 2]]]
    n.g2 <- as.character(names(g.2))
    #n.2 <- length(g.2)
    if(!shared) {
      n.g1 <- n.g1[!n.g1 %in% n.g2]
      n.g2 <- n.g2[!n.g2 %in% n.g1]
    }
    
    tab <- gene_network_adj[n.g1, n.g2, drop = F]
    out <- list(tab, n.g1, n.g2)
    
    return(out)
  }, mc.cores = mc_cores_ct)
  
  
  idx <- which(mclapply(sub_adj, function(n) sum(n[[1]], na.rm = T), mc.cores = mc_cores_ct) !=0)
  
  if(length(idx) == 0 ) {
    print("No CT available")
    return("No CT available")
  } else {
    
    sub_adj <- sub_adj[idx]
    comb_p <- comb_p[idx,]
    perm_list <- mclapply(1:nrow(comb_p), function(x) {
      n.1 <- l_cl_list[comb_p[x,1]]
      n.2 <- l_cl_list[comb_p[x,2]]
      out <- perm_link(r = n.1, c = n.2, gene_network_adj, core = mc_cores_perm, k, hash = F)
      return(out)
    }, mc.cores = mc_cores_ct)
    
    len1 <- array(unlist(mclapply(sub_adj, function(n) length(n[[2]]), mc.cores = mc_cores_ct)), 
                  dim = c(length(sub_adj), 1))
    len2 <- array(unlist(mclapply(sub_adj, function(n) length(n[[3]]), mc.cores = mc_cores_ct)), 
                  dim = c(length(sub_adj), 1))
    comb_p_len <- matrix(cbind(len1, len2), ncol = 2)
    unq_comb_p_len <- unique(comb_p_len)
    perm_list <- mclapply(1:nrow(unq_comb_p_len), function(x) {
      n.1 <- unq_comb_p_len[x, 1]
      n.2 <- unq_comb_p_len[x, 2]
      out <- perm_link(r = n.1, c = n.2, gene_network_adj, core = mc_cores_perm, k, hash = hash)
      return(out)
    }, mc.cores = mc_cores_ct)
    
    unq_comb_p_len <- data.frame(unq_comb_p_len)
    unq_comb_p_len$var <- 1:nrow(unq_comb_p_len)
    unq_comb_p_len <- acast(unq_comb_p_len, X1 ~ X2, value.var = "var")
    
    
    ans <- parallel::mclapply(1:nrow(comb_p), function(z) {
      tab <- sub_adj[[z]][[1]]
      n.1 <- length(sub_adj[[z]][[2]])
      n.2 <- length(sub_adj[[z]][[3]])
      nlink <- matrix(sum(tab))
      idx <- unq_comb_p_len[as.character(n.1), as.character(n.2)]
      perm_l <- c(list(nlink), perm_list[[idx]])
      p_val <- calc_p(perm_l)
      
      g.1 <- gs_list[[comb_p[z, 1]]]
      g.2 <- gs_list[[comb_p[z, 2]]]
    
      ct <- cross_talk(mat = tab, weight = list(g1 = g.1, g2 =g.2))
      out <- array(c(comb_p[z, 1], comb_p[z, 2], ct[1:4], as.vector(p_val), ct[5:8]), dim = c(1, 11))
      
      return(out)
    }, mc.cores = mc_cores_perm)
    
    ans <- do.call(rbind, ans)
    all.v <- unlist(perm_list)
    all.v <- c(as.numeric(ans[,6]), all.v)
    link_FDR <- eFDR(real_values = as.numeric(ans[,6]), all_values = all.v, mc.cores = mc_cores_ct)
    p_val_BH <- stats::p.adjust(as.numeric(ans[,7]), method = "BH")
    ans <- data.frame(gs1 = ans[, 1], 
                      gs2 = ans[, 2], 
                      ct_score = as.numeric(ans[, 3]), 
                      ngenes_gs1 = as.numeric(ans[, 4]), 
                      ngenes_gs2 = as.numeric(ans[, 5]),
                      nlink = as.numeric(ans[, 6]), 
                      p_value_link = as.numeric(ans[, 7]), 
                      FDR_link = link_FDR, 
                      p_adj_BH = p_val_BH, 
                      weight_gs1 = as.numeric(ans[, 8]), 
                      weight_gs2 = as.numeric(ans[, 9]), 
                      genes_gs1 = ans[, 10], 
                      genes_gs2 = ans[, 11], 
                      stringsAsFactors = F)
    if(ct_info) {
      ct_info <- parallel::mclapply(1:nrow(comb_p), function(z) {
        cl1 <- gs_list[[comb_p[z, 1]]]
        ncl1 <- names(cl1)
        cl2 <- gs_list[[comb_p[z, 2]]]
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
      },mc.cores = mc_cores_ct)
      ct_info <- do.call(rbind, ct_info)
      colnames(ct_info) <- c("gs1", "gs1_gene", "gs2", "gs2_gene", "score")
      ct_info <- data.frame(gs1 = ct_info[, 1],
                            gs1_gene = ct_info[, 2],
                            gs2 = ct_info[, 3],
                            gs2_gene = ct_info[, 4],
                            score = as.numeric(ct_info[, 5]),
                            stringsAsFactors = F)
      
      return(list(ct_info = ct_info, ct_res= ans))
    } else {
      return(ans)
    }
    
  }
  
  
}

