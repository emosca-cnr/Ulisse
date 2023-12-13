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

plot_gsea_heatmap <- function(gsea_res=NULL, nes_sign=TRUE, a=0.25, p.stat="FDRq", na_col="khaki", max_gs=50, min.p=NA, ...){
  
  if(length(gsea_res)<2){
    stop("This function requires at least two GSEAs.\n")
  }
  
  if(!p.stat %in% colnames(gsea_res[[1]])){
    stop("p.stat", p.stat, "is not in", colnames(ora_res[[1]]), ".\n")
  }
  
  cat("Size: ", unlist(lapply(gsea_res, nrow)), "\n")
  
  X_matrix <- merge(gsea_res[[1]][, c("description", p.stat)], gsea_res[[2]][, c("description", p.stat)], by=1, all=T, sort=F)
  colnames(X_matrix)[2:3] <- names(gsea_res)[1:2]
  
  if(length(gsea_res)>2){
    for(i in 3:length(gsea_res)){
      X_matrix <- merge(X_matrix, gsea_res[[i]][, c("description", p.stat)], by=1, all = T, sort=F)
      colnames(X_matrix)[i+1] <- names(gsea_res)[i]
    }
  }
  
  rownames(X_matrix) <- X_matrix$description
  X_matrix$description <- NULL
  X_matrix <- -log10(X_matrix)

  if(is.null(na_col)){
    X_matrix[is.na(X_matrix)] <- 0
  }else{
    X_matrix[X_matrix < -log10(a)] <- NA
  }
  X_matrix <- as.matrix(X_matrix)
  
  if(nrow(X_matrix)>max_gs){
    cat("Only the top", max_gs, "gene sets will be considered...\n")
    X_matrix <- X_matrix[1:max_gs, ]
  }
  
  if(any(is.na(X_matrix))){
    cat("The presence of many NA values may cause errors in hclust(). To avoid this error, (i) increase 'a' or (ii) set cluster_columns = F and/or cluster_rows = F.\n")
  }
  
  if(nes_sign){
    
    nes_sign <- merge(gsea_res[[1]][, c("description", "nes")], gsea_res[[2]][, c("description", "nes")], by=1, all = T, sort=F)
    colnames(nes_sign)[2:3] <- names(gsea_res)[1:2]
    
    if(length(gsea_res)>2){
      for(i in 3:length(gsea_res)){
        nes_sign <- merge(nes_sign, gsea_res[[i]][, c("description", "nes")], by=1, all = T, sort=F)
        colnames(nes_sign)[1+i] <- names(gsea_res)[i]
      }
    }
    rownames(nes_sign) <- nes_sign$description
    nes_sign$description <- NULL
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
  
  if(is.na(min.p)){
    min.p <- min(10^-(X_matrix), na.rm = T)
  }else{
    X_matrix[X_matrix > -log10(min.p)] <- -log10(min.p)
  }
  
  col_fun <- colorRamp2(seq(min(X_matrix, na.rm = T), min(max(X_matrix, na.rm=T), -log10(min.p)), length.out=5), brewer.purples(5))
  
  plot(Heatmap(X_matrix, col=col_fun, na_col = na_col, name=paste0("-log10(", p.stat, ")"), rect_gp = gpar(col = "black", lwd = 1), cell_fun = function(j, i, x, y, width, height, fill) { grid.text(nes_sign[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold", col=ifelse(nes_sign[i, j]=="+", "limegreen", "firebrick")))}, ...))
  
}