#' Function to caluclate sub-components pathway score
#' @description Calculates the score of each pathway connected components
#' @details By coupling the pathway data to the adjacency matrix, the function identifies the pathway components as
#'  connected components (see [igraph::components()] for details). For each component, the function calcualtes the score 
#'  as described in the article.
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj adjacency matrix of the whole gene network considered
#' @param weight weights of the genes in pathway list. If not provided, the function assign to each gene
#' a weight of 1
#' @param mc_cores_cc numebr of threads to be used for cc calculation
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
#' }
#' }
#' @examples 
#' \dontrun{
#' ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#' adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#' ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#' wgt <- rep(1, 6)
#' p_cc <- pathway_cc(pathway_list = ptw_list, 
#' gene_network_adj = adj, weight = wgt, mc_cores_cc = 1)
#' }
#' @import parallel
#' @import igraph
#' @importFrom stringi stri_c
#' @export


pathway_cc <- function (pathway_list, gene_network_adj,
                        weight, mc_cores_cc = 2) {
  genes <- rownames(gene_network_adj)
  if(is.null(weight) ) {
    weight <- rep(1, length(genes))
    names(weight) <- genes
  } 
  if(!is(gene_network_adj, "sparseMatrix" )) {
    gene_network_adj <- as(gene_network_adj, "dgCMatrix")
  }
  gene_network_adj <- sign(gene_network_adj)
  weight <- weight[weight !=0]
  sub_adj_mt <- gene_network_adj[as.character(names(weight)), as.character(names(weight)), drop = F]
  
  
  xx <- mclapply(pathway_list, function(x) {
    x <- as.character(x)
    tmp <- sub_adj_mt[x, x, drop = F]
    
    return(tmp)
  }, mc.cores = mc_cores_cc)
  
  
  idx <- which(lapply(xx, function(n) sum(n)) !=0)
  if(length(idx) == 0) {
    print("No pathway on which calculate CC")
    return("No pathway on which calculate CC")
  } else {
    xx <- xx[idx]
    pathway_list <- pathway_list[idx]
    
    
    cc_out <- mclapply(1:length(xx), function(x) {
      tab <- xx[[x]]
      cc <- igraph::components(graph = igraph::graph_from_adjacency_matrix(tab,
                                                                           mode = "undirected"))
      
      
    },   mc.cores = mc_cores_cc)
    
    out <- mclapply(1:length(pathway_list), function(z) {
      tab <- xx[[z]]
      name.tab <- names(xx)[z]
      cc <- cc_out[[z]]
      tmp <- lapply(1:cc$no, function(x){
        cc_genes <- as.character(names(cc$membership[which(cc$membership==x)]))
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
    }, mc.cores = mc_cores_cc)
    out <- do.call(rbind, out)
    idx <- out[, "score"]== 0
    out <- out[!idx,]
    
    cc_list <- lapply(cc_out, function(x) {
      x <- x$membership
      return(x)
    })
    
    return(list(membership = cc_list, 
                pathway_cc = out))
  }
  
}

