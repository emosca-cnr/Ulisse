#' Function to calculate sub-components pathway score
#' @description The function calculates the score of each pathway connected components
#' @details By coupling the pathway data to the adjacency matrix, the function identifies the pathway components as
#' connected components (see [igraph::components()] for details). For each component, the function calculates the score 
#' as described in the article.
#' @param gs_list a named list of genes grouped into gene-sets as obtained from `preparing_gs_list`, `preparing_msigdb_list()`, 
#'  `preparing_DEG_list()` or `preparing_expr_list()` functions
#' @param gene_network_adj adjacency matrix of the whole gene network considered
#' @param mc_cores_cc number of threads to be used for cc calculation
#' @return The function returns a list of two object:
#' \enumerate{
#'  \item membership: results of [igraph::components()] for each pathway considered
#'  \item {pathway_cc}: a data frame with
#' \itemize{
#'  \item pathway: name of the pathway where the components scores are calculated
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
#' @importFrom methods is as
#' @export


gene_set_cc <- function (gs_list, gene_network_adj, mc_cores_cc = 2) {
  
  if(!is(gene_network_adj, "sparseMatrix" )) {
    gene_network_adj <- as(gene_network_adj, "dgCMatrix")
  }
  gene_network_adj <- sign(gene_network_adj)
  gs_list <- lapply(gs_list, function(x) x <- x[names(x) %in% rownames(gene_network_adj)])
  gene <- as.character(unique(unlist(lapply(gs_list, names))))
  gene_network_adj <- gene_network_adj[gene, gene, drop = F]
  
  xx <- mclapply(gs_list, function(x) {
    
    tmp <- gene_network_adj[as.character(names(x)), as.character(names(x)), drop = F]
    
    return(tmp)
  }, mc.cores = mc_cores_cc)
  
  
  idx <- which(lapply(xx, function(n) sum(n)) !=0)
  if(length(idx) == 0) {
    print("No pathway on which calculate CC")
    return("No pathway on which calculate CC")
  } else {
    xx <- xx[idx]
    gs_list <- gs_list[idx]
    
    
    cc_out <- mclapply(1:length(xx), function(x) {
      tab <- xx[[x]]
      cc <- igraph::components(graph = igraph::graph_from_adjacency_matrix(tab,
                                                                           mode = "undirected"))
      return(cc)
      
    },   mc.cores = mc_cores_cc)
    
    out <- mclapply(1:length(gs_list), function(z) {
      tab <- xx[[z]]
      name.tab <- names(xx)[z]
      cc <- cc_out[[z]]
      tmp <- lapply(1:cc$no, function(x){
        cc_genes <- as.character(names(cc$membership[which(cc$membership==x)]))
        tab.tmp <- tab[cc_genes, cc_genes, drop = F]
        score <- sum(t(gs_list[[z]][cc_genes]) %*% tab.tmp)/2
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
    out <- data.frame(pathway = out[,1],
                      ID = out[, 2],
                      score = as.numeric(out[,3]),
                      n_gene = as.numeric(out[,4]),
                      n_link = as.numeric(out[,5]),
                      gene = out[,5], stringsAsFactors = F)
    
    cc_list <- lapply(cc_out, function(x) {
      x <- x$membership
      return(x)
    })
    
    return(list(membership = cc_list, 
                pathway_cc = out))
  }
  
}

