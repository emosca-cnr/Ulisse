#' Plot heatmap of ORA results for multiple runs
#' @param ora_res result of function gsea()
#' @param nes_sign whether to print nes sign
#' @param a alpha over FDR
#' @param na_col color for FDR > a
#' @param max_gs maximum number of gene sets that will be plotted
#' @param ... further arguments to ComplexHeatmap::Heatmap()
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom pals brewer.purples
#' @importFrom grid gpar grid.text

plot_ora_heatmap <- function(ora_res=NULL, p.stat="p_adj", a=0.25, na_col="khaki", max_gs=50, ...){
  
  if(length(ora_res)<2){
    stop("This function requires at least two ORAs.\n")
  }
  
  cat("Size: ", unlist(lapply(ora_res, nrow)), "\n")
  
  if(!p.stat %in% colnames(ora_res[[1]])){
    stop("p.stat", p.stat, "is not in", colnames(ora_res[[1]]), ".\n")
  }
  
  X_matrix <- merge(ora_res[[1]][, c("id", p.stat)], ora_res[[2]][, c("id", p.stat)], by=1, all=T, sort=F)
  colnames(X_matrix)[2:3] <- names(ora_res)[1:2]
  
  if(length(ora_res)>2){
    for(i in 3:length(ora_res)){
      X_matrix <- merge(X_matrix, ora_res[[i]][, c("id", p.stat)], by=1, all = T, sort=F)
      colnames(X_matrix)[i+1] <- names(ora_res)[i]
    }
  }
  
  rownames(X_matrix) <- X_matrix$id
  X_matrix$id <- NULL
  #rownames(X_matrix) <- gss$gs_name[match(rownames(X_matrix), gss$gs_id)]
  X_matrix <- -log10(X_matrix)
  X_matrix[X_matrix < -log10(a)] <- NA
  X_matrix <- as.matrix(X_matrix)
  
  if(nrow(X_matrix)>max_gs){
    cat("Only the top", max_gs, "gene sets will be considered...\n")
    X_matrix <- X_matrix[1:max_gs, ]
  }
  
  if(any(is.na(X_matrix))){
    warning("The presence of many NA values may cause errors in hclust(). To avoid this error, (i) increase 'a' or (ii) set cluster_columns = F and/or cluster_rows = F.\n")
  }
  
  Heatmap(X_matrix, col=brewer.purples(5), na_col = na_col, name="-log10(p)", rect_gp = gpar(col = "black", lwd = 1), ...)
  
}
