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
#' @export


pathway_cc <- function (pathway_list, gene_network_adj,
                        weight, k=9, mc_cores_cc = 2, mc_cores_perm = 2
) {
  gene_network_adj <- sign(gene_network_adj)
  perm_list <- mclapply(1:k, function(x) {
    tmp <- matrix(as.numeric(gene_network_adj), ncol = ncol(gene_network_adj), 
                  dimnames = list(sample(rownames(gene_network_adj), nrow(gene_network_adj))))
    colnames(tmp) <- rownames(tmp)
    return(tmp)
  }, mc.cores = mc_cores_perm)
  
  xx <- mclapply(pathway_list, function(x) {
    idx <- which(rownames(gene_network_adj) %in% x)
    return(gene_network_adj[idx,idx, drop = F])
  }, mc.cores = mc_cores_cc)
  names(weight) <- rownames(gene_network_adj)
  
  
  idx <- which(lapply(xx, function(n) sum(n)) !=0)
  if(length(idx) == 0) {
    print("No pathway on which calculate CC")
    return("No pathway on which calculate CC")
  } else {
    xx2 <- xx[idx]
    
    xxPL <- mclapply(1:k, function(n){
      mclapply(pathway_list[idx], function(x) {
        idx <- which(rownames(perm_list[[n]]) %in% x)
        return(perm_list[[n]][idx,idx, drop = F])
      }, mc.cores = mc_cores_cc)
    }, mc.cores = mc_cores_perm )
    
    
    pct <- mclapply(1:length(xx2), function(x) {
      
      cc <- igraph::components(graph = igraph::graph_from_adjacency_matrix(as.matrix(xx2[[x]]),
                                                                           mode = "undirected"))
      cc_data <- data.frame(pathway = rep(names(xx2)[x], cc$no),
                            ID = 1:cc$no,
                            score = rep(0, cc$no),
                            n_gene = rep(0, cc$no),
                            n_link = rep(0, cc$no), 
                            gene = rep(0, cc$no),
                            stringsAsFactors = F)
      
      
      for(j in 1:cc$no) {
        cc_genes <- names(cc$membership[which(cc$membership==j)])
        gene_network_adj_ijW <- sum(t(weight[cc_genes]) %*% as.matrix(xx2[[x]][cc_genes, cc_genes]))/2
        gene_network_adj_ijW <- gene_network_adj_ijW/length(cc_genes)
        cc_data$score[j] <- gene_network_adj_ijW
        cc_data$n_gene[j] <- length(cc_genes)
        cc_data$n_link[j] <- sum(xx2[[x]][cc_genes, cc_genes])/2
        cc_data$gene[j] <- paste(cc_genes, collapse = ";")
        
        
      }
      perm_cc <- mclapply(1:length(xxPL), function(j) {
        tmp <- matrix(data=0, nrow = nrow(cc_data), dimnames = list(cc_data$ID))
        z = 1
        for (k in cc_data$ID) {
          cc_genes <- names(cc$membership[which(cc$membership==k)])
          gene_network_adj_ijW <- sum(t(weight[cc_genes]) %*% as.matrix(xxPL[[j]][[x]][cc_genes, cc_genes]))
          gene_network_adj_ijW <- gene_network_adj_ijW/length(cc_genes)
          tmp[z,1] <- gene_network_adj_ijW
          z = z+1
        }
        return(tmp)
      }, mc.cores = mc_cores_perm)
      
      tmp <- cc_data[,3, drop=F]
      rownames(tmp) <- cc_data$ID
      perm_cc <- c(list(tmp), perm_cc)
      cc_data$p_value <- calc_p(perm_cc)
      
      cc_list <- list(components_results = cc,
                      pathway_cc = cc_data,
                      perm = perm_cc)
      zero <- cc_list$pathway_cc$ID[which(cc_list$pathway_cc$score == 0)]
      if(length(zero) != 0) {
        cc_list$components_results$membership <- cc_list$components_results$membership[-which(cc_list$components_results$membership %in% zero)]
        cc_list$pathway_cc <- cc_list$pathway_cc[which(cc_list$pathway_cc$score >0),]
        cc_list$perm <- lapply(cc_list$perm, function(l) {
          l <- l[-zero]
          return(l)
        })
      }
      
      return(cc_list)
      
    },   mc.cores = mc_cores_cc)
    names(pct) <- names(xx2)
    
    
    
    
    cc_list <- mclapply(1:length(pct), function (k) {
      tmp <- pct[[k]]$components_results$membership
      return(tmp)
    }, mc.cores = mc_cores_cc)
    names(cc_list) <- names(xx2)
    
    cc_data <- mclapply(1:length(pct), function(k) {
      tmp <- pct[[k]]$pathway_cc
      return(tmp)
    }, mc.cores = mc_cores_cc)
    cc_data <- do.call(rbind, cc_data)
    
    perm <- mclapply(1:length(pct), function(k) {
      tmp <- pct[[k]]$perm
      return(tmp)
    }, mc.cores = mc_cores_cc)
    
    cc_data$eFDR <- eFDR(real_values = as.vector(unlist(cc_data$score)), 
                         all_values = as.vector(unlist(perm)))
    
    return(list(membership = cc_list, 
                pathway_cc = cc_data))
  }
  
}

