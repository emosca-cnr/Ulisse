#' Plot heatmap of GSEA results for multiple runs
#' @param gsea_res result of function gsea()
#' @param nes_sign whether to print nes sign
#' @param a alpha over FDR
#' @param na_col color for FDR > a
#' @param max_gs maximum number of gene sets that will be plotted
#' @param ... further arguments to ComplexHeatmap::Heatmap()
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom pals brewer.purples
#' @importFrom grid gpar grid.text

plot_gsea_heatmap <- function(gsea_res=NULL, nes_sign=TRUE, a=0.25, na_col="khaki", max_gs=50, ...){
  
  if(length(gsea_res)<2){
    stop("This function requires at least two GSEAs.\n")
  }
  
  cat("Size: ", unlist(lapply(gsea_res, nrow)), "\n")
  
  X_matrix <- merge(gsea_res[[1]][, c("id", "FDRq")], gsea_res[[2]][, c("id", "FDRq")], by=1, all=T, sort=F)
  colnames(X_matrix)[2:3] <- names(gsea_res)[1:2]
  
  if(length(gsea_res)>2){
    for(i in 3:length(gsea_res)){
      X_matrix <- merge(X_matrix, gsea_res[[i]][, c("id", "FDRq")], by=1, all = T, sort=F)
      colnames(X_matrix)[i+1] <- names(gsea_res)[i]
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
  
  if(nes_sign){
    
    nes_sign <- merge(gsea_res[[1]][, c("id", "nes")], gsea_res[[2]][, c("id", "nes")], by=1, all = T, sort=F)
    colnames(nes_sign)[2:3] <- names(gsea_res)[1:2]
    
    if(length(gsea_res)>2){
      for(i in 3:length(gsea_res)){
        nes_sign <- merge(nes_sign, gsea_res[[i]][, c("id", "nes")], by=1, all = T, sort=F)
        colnames(nes_sign)[1+i] <- names(gsea_res)[i]
      }
    }
    rownames(nes_sign) <- nes_sign$id
    nes_sign$id <- NULL
    nes_sign <- sign(nes_sign)
    nes_sign[which(nes_sign>0, arr.ind = T)] <- "+"
    nes_sign[which(nes_sign==-1, arr.ind = T)] <- "-"
    nes_sign[which(is.na(nes_sign), arr.ind = T)] <- ""
  }else{
    nes_sign <- X_matrix
    nes_sign[] <- ""
  }
  
  nes_sign <- nes_sign[match(rownames(X_matrix), rownames(nes_sign)), ]
  nes_sign[is.na(X_matrix)] <- ""
  
  Heatmap(X_matrix, col=brewer.purples(5), na_col = na_col, name="-log10(q)", rect_gp = gpar(col = "black", lwd = 1), cell_fun = function(j, i, x, y, width, height, fill) { grid.text(nes_sign[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold", col=ifelse(nes_sign[i, j]=="+", "firebrick", "limegreen")))}, ...)
  
}