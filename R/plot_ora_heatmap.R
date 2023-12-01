#' Plot heatmap of ORA results for multiple runs
#' @param ora_res result of function gsea()
#' @param a alpha over FDR
#' @param na_col color for FDR > a
#' @param max_gs maximum number of gene sets that will be plotted
#' @param p.stat p-value that will be plotted
#' @param ... further arguments to ComplexHeatmap::Heatmap()
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom pals brewer.purples
#' @importFrom grid gpar grid.text
#' @importFrom circlize colorRamp2

plot_ora_heatmap <- function(ora_res=NULL, p.stat="p_adj", a=0.25, na_col="khaki", min.p=0.0001, max_gs=50, ...){
  
  if(length(ora_res)<2){
    stop("This function requires at least two ORAs.\n")
  }
  
  cat("Size: ", unlist(lapply(ora_res, nrow)), "\n")
  
  if(!p.stat %in% colnames(ora_res[[1]])){
    stop("p.stat", p.stat, "is not in", colnames(ora_res[[1]]), ".\n")
  }
  
  X_matrix <- merge(ora_res[[1]][, c("description", p.stat)], ora_res[[2]][, c("description", p.stat)], by=1, all=T, sort=F)
  colnames(X_matrix)[2:3] <- names(ora_res)[1:2]
  
  if(length(ora_res)>2){
    for(i in 3:length(ora_res)){
      X_matrix <- merge(X_matrix, ora_res[[i]][, c("description", p.stat)], by=1, all = T, sort=F)
      colnames(X_matrix)[i+1] <- names(ora_res)[i]
    }
  }
  
  rownames(X_matrix) <- X_matrix$description
  X_matrix$description <- NULL
  #rownames(X_matrix) <- gss$gs_name[match(rownames(X_matrix), gss$gs_id)]
  X_matrix <- -log10(X_matrix)
  
  if(is.null(na_col)){
    X_matrix[is.na(X_matrix)] <- 0
  }else{
    X_matrix[X_matrix < -log10(a)] <- NA
  }
  X_matrix <- as.matrix(X_matrix)
  
  if(nrow(X_matrix)>max_gs){
    cat("Found", nrow(X_matrix), "gene sets, only the top", max_gs, "gene sets will be considered...\n")
    X_matrix <- X_matrix[1:max_gs, ]
  }
  
  if(any(is.na(X_matrix))){
    warning("The presence of many NA values may cause errors in hclust(). To avoid this error, (i) increase 'a' or (ii) set cluster_columns = F and/or cluster_rows = F.\n")
  }
  col_fun <- colorRamp2(seq(0, -log10(min.p), length.out=5), brewer.purples(5))
  plot(Heatmap(X_matrix, col=col_fun, na_col = na_col, name=paste0("-log10(", p.stat, ")"), rect_gp = gpar(col = "black", lwd = 1), ...))
  
}
