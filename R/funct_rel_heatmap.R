#' Function to obtain a detailed visualization of specific results of functional relevance analysis
#' @description The functions highlight the links between defined genes and/or pathways to better inspect the results obtained with `gene_functional_relevance`
#' analysis
#' @details The function uses the same input used for `gene_functional_relevance` to retrieve the gene interactions of specific `target_g` and/or `target_ptw`.
#' The function produce an heatmap with the genes that interacts with the user-provided ones for each pathway considered. In the heatmap can be reported the 
#' product of the weight of the genes if `weight` is provided, or the existence of the link (0-1) otherwise.
#' @param pct output of the `pathway_cross_talk()` function
#' @param adj adjacency matrix used for PCT calculation
#' @param target_g vector with the genes of interest. If provided, the function will consider only the gene that interact with these in the plot
#' @param target_ptw vector with the pathways of interest. If provided, the function will consider only the CT where these pathways are involved
#' @param weight weight of the genes. If provided, the function will plot the product of the weights in the heatmpa. If `NULL` the function will lot only
#' the presence/absence of the links (1-0)
#' @param n_cores number of cores to be used by the function
#' @param colors vector of two colors to be used in the heatmap. If not provided, the functions used ble and red
#' @param width,height,res,units graphical value of `jpeg()` function
#' @param ... further arguments passed to `Heatmap` function
#' @importFrom grDevices dev.off jpeg
#' @importFrom circlize colorRamp2
#' @import igraph
#' @import complexHeatmap
#' @export


funct_rel_heatmap <- function(pct, adj, target_g = NULL, target_ptw = NULL,
                              weight = wgth, n_cores = 1, file_name = NULL, 
                              colors = NULL, width = 200, height = 200, res = 300, units = "mm", ...) {
  
  if(!is.null(target_g)) {
    idx1 <- grep(paste(target_g, collapse =  "|"), pct$gene_pathway1)
    idx2 <- grep(paste(target_g, collapse =  "|"), pct$gene_pathway2)
    idx <- union(idx1, idx2)
    pct <- pct[idx,]
  }
  
  if(!is.null(target_ptw)) {
    idx1 <- grep(paste(target_ptw, collapse =  "|"), pct$pathway_1)
    idx2 <- grep(paste(target_ptw, collapse =  "|"), pct$pathway_2)
    idx <- union(idx1, idx2)
    pct <- pct[idx,]
  }
  
  if(is.null(weight) ) {
    weight <- rep(1, length(genes))
    names(weight) <- genes
  } 
  
  if(!is(adj, "sparseMatrix" )) {
    adj <- as(adj, "dgCMatrix")
  }
  out <- mclapply(1:nrow(pct), function(x) {
    tmp <- pct[x,]
    g1 <- unlist(strsplit(tmp$gene_pathway1, ";", fixed = T))
    g2 <- unlist(strsplit(tmp$gene_pathway2, ";", fixed = T))
    sub.adj <- as.matrix(adj[g1, g2, drop = F])
    sub.adj <- reshape2::melt(sub.adj)
    sub.adj$w <- weight[sub.adj$Var1] * weight[sub.adj$Var2]
    sub.adj$Var1 <- paste(rep(tmp$pathway_1, nrow(sub.adj)), sub.adj$Var1, sep = "_")
    sub.adj$Var2 <- paste(rep(tmp$pathway_2, nrow(sub.adj)), sub.adj$Var2, sep = "_")
    return(sub.adj)
  }, mc.cores = n_cores)
  out <- do.call(rbind, out)
  out.g <- graph_from_data_frame(out, directed = F)
  out2 <- as_adjacency_matrix(out.g, sparse = F, attr = "w")
  
  
  if(is.null(file_name)) {
    file_name <- "funct_rel_heatmap.jpeg"
  }
  
  if(is.null(colors)) {
    colors <- c("lightyellow", "red")
  }
  colors <-  colorRamp2(c(min(out2), max(out2)), colors)
  jpeg(filename = file_name, width = width, height = height, res = res, units = units)
  h <- Heatmap(as.matrix(out2), heatmap_legend_param = list(title="Score"),
               col = colors, ...)
  draw(h, heatmap_legend_side = "left")
  dev.off()
  return(out2)
}
