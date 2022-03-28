#' Function to caluclate sub-components pathway score
#' @description Calculates the score of each pathway connected components
#' @details By coupling the pathway data to the adjacency matrix, the function identifies the pathway components as
#'  connected components (see [igraph::components()] for details). For each component, the function calcualtes the score 
#'  as described in the article. Ths approach is applied also to all permuted version of the adjacency matrix. 
#'  The score and the permuted score are then used to calcualte p-value and FDR
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param weight weights of the genes in gene_network_adj
#' @param mc_cores_cc numebr of threads to be used for cc calculation
#' @param mc_cores_perm number of threads to be used for permutations
#' @param k numebr of permutation of the adjacency matrix
#' @return The function returns a list of two object:
#' \enumerate{
#'  \item membership: results of [igraph::components()] for each pathway considered
#'  \item {pathway_cc}: a data frame with
#' \itemize{
#'  \item pathway: name of the pathway where the components scores are calcualted
#'  \item ID: id of the component of the pathway
#'  \item score: component score
#'  \item n_gene: number of genes in the component
#'  \item n_link: number of links in the component
#'  \item gene: names of the genes in the component
#'  \item p_value: empirical p-value calculated by using the permutation approcah
#'  \item eFDR: empirical FDR
#' }
#' }
#' @examples 
#' ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#' adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#' ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#' wgt <- rep(1, 6)
#' p_cc <- pathway_cc(pathway_list = ptw_list, 
#' gene_network_adj <- adj, weight = wgt, mc_cores_cc = 1, mc_cores_perm = 1, k = 9)
#' @import parallel
#' @import igraph
#' @import stringi
#' @export


pathway_cc <- function (pathway_list, gene_network_adj,
                        weight, k=9, mc_cores_cc = 2, mc_cores_perm = 2
) {
  gene_network_adj <- sign(gene_network_adj)
  names(weight) <- rownames(gene_network_adj)
  weight <- weight[weight !=0]
  gene_network_adj[names(weight), names(weight)]
  
  
  xx <- mclapply(pathway_list, function(x) {
    tmp <- gene_network_adj[x, x, drop = F]
    
    return(tmp)
  }, mc.cores = mc_cores_cc)
  
  
  idx <- which(lapply(xx, function(n) sum(n)) !=0)
  if(length(idx) == 0) {
    print("No pathway on which calculate CC")
    return("No pathway on which calculate CC")
  } else {
    xx <- xx[idx]
    pathway_list <- pathway_list[idx]
    
    perm_list <- mclapply(1:k, function(x) {
      tmp <- gene_network_adj
      rownames(tmp) <- sample(rownames(tmp), nrow(tmp))
      colnames(tmp) <- rownames(tmp)
      out <- mclapply(pathway_list, function(x) {
        tmp.2 <- tmp[x, x, drop = F]
        
        return(tmp.2)
      }, mc.cores = mc_cores_cc)
    }, mc.cores = mc_cores_perm)
    
    
    
    cc_out <- mclapply(1:length(xx), function(x) {
      tab <- xx[[x]]
      cc <- igraph::components(graph = igraph::graph_from_adjacency_matrix(tab,
                                                                           mode = "undirected"))
      
      
    },   mc.cores = mc_cores_cc)
    
    all <- c(list(xx), perm_list)
    
    out <- mclapply(1:(k+1), function(j) {
      out.1 <- mclapply(1:length(pathway_list), function(z) {
        tab <- all[[j]][[z]]
        name.tab <- names(all[[j]])[z]
        cc <- cc_out[[z]]
        tmp <- lapply(1:cc$no, function(x){
          cc_genes <- names(cc$membership[which(cc$membership==x)])
          tab.tmp <- tab[cc_genes, cc_genes, drop = F]
          score <- sum(t(weight[cc_genes]) %*% tab.tmp)/2
          score <- score/length(cc_genes)
          n_gene <- length(cc_genes)
          n_link <- sum(tab.tmp)/2
          gene.name <- stri_c(cc_genes, collapse = ";")
          tab.out <- matrix(c(name.tab, x, score, n_gene, n_link, gene.name), nrow = 1)
        })
        tmp <- do.call(rbind, tmp)
        colnames(tmp) <- c("pathway", "ID", "score", "n_gene", "n_link", "gene")
        return(tmp)
      })
      out.1 <- do.call(rbind, out.1)
    })
    idx <- out[[1]][, "score"]== 0
    out <- lapply(out, function(x) {
      x <- x[!idx,]
      return(x)
    })
    p_list <- mclapply(out, function(x) {
      tmp <- matrix(as.numeric(x[,"score"]), ncol = 1)
      rownames(tmp) <- x[, 1]
      return(tmp)
    }, mc.cores = mc_cores_perm)
    p_val <- calc_p(p_list)
    
    
    cc_list <- lapply(cc_out, function(x) {
      x <- x$membership
      return(x)
    })
    
    eFDR <- eFDR(real_values = as.vector(unlist(p_list[[1]])), 
                 all_values = as.vector(unlist(p_list)))
    out <- data.frame(out[[1]], stringsAsFactors = F)
    out$p_value <- p_val[,1]
    out$FDR <- eFDR
    return(list(membership = cc_list, 
                pathway_cc = out))
  }
  
}

