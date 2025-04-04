#' Cell-cell communication network plot
#' @description This function is used to plot the results of the cell-cell communication analysis
#'  overlapped to the single-cell plot
#' @details This functions uses the cell-embeddings produced by the single-cell analysis
#'  and the cell-cell communication analysis results, obtained with `gs_cross_talk()`, to plot
#'  the communication network overlapped to the cell plot. This functions do not filter the 
#'  `gs_cross_talk()` results to plot only the significant ones.
#'  Cells and communication network vertices are colored according to the annotation, provided with `cl_cell`
#' @param cell_emb cell embedding (UMAP or tSNE) matrix
#' @param cl_res cell-cell communication results (from `gs_cross_talk()`) 
#' @param cl_cell two column data-frame, the first should be the cell names, the second the annotation (cluster or cell types)
#' @param palette named vector of the color associated to each cluster or cell types, and it will be used to set both 
#'  cell points and communication network vertices color. If NULL ignored, ggplot automatic palette used
#' @param e_scale value to scale the edge weights for plotting, change to adjust maximum link widths
#' @param edge_alpha = 0.8 transparency of the edge color
#' @param edge_color = edge color
#' @param label logical, if cluster/cell type name should be plotted on the vertex names
#' @param save_file saving plot name. If NULL ignored
#' @param point_size = 0.5 size of the cell points to be plotted.
#' @param point_alpha = 0.5 transparency of the color fill of the cell points
#' @param node_size = 5 size of the vertices of the cell-cell communication network overlapped to the cell points
#' @param text_size size of the node labels plotted (if `label = T`)
#' @param text_color = text color (if `label = T`)
#' @param heigh,width,res,unit params used to save plot (in jpeg format)
#' @importFrom grDevices dev.off jpeg boxplot.stats
#' @importFrom igraph graph_from_data_frame get.edge.ids
#' @importFrom ggplot2 ggplot aes geom_point guides theme_test  scale_color_manual scale_fill_manual
#' @importFrom ggraph geom_edge_link geom_node_point geom_node_text
#' @export


plot_communication_network <- function(cell_emb, cl_res, cl_cell, palette, e_scale = 9, label, 
                                       save_file, point_size = 0.5, point_alpha = 0.5, 
                                       edge_alpha = 0.8, edge_color = "gray45",
                                       node_size = 5, text_size = 2, text_color = "black", heigh = 100, width = 100, res = 300, unit = "mm") {
  
  colnames(cl_cell) <- c("cell", "cluster")
  
  cell_emb <- data.frame(cell_emb)
  cell_emb <- cell_emb[cl_cell$cell,]
  cell_emb <- cell_emb[, 1:2]
  #colnames(cell_emb) <- c("dim_1", "dim_2")
  cell_emb <- cbind(cl_cell, cell_emb)
  n <- unique(cell_emb$cluster)
  
  cl_coord <- lapply(n, function(x) {
    tmp <- cell_emb[cell_emb$cluster == x,]
    
    tmp.x <- tmp[, 3]
    s <- boxplot.stats(tmp.x)$stats
    tmp[, 3][tmp.x < s[1]] <- s[1]
    tmp[, 3][tmp.x > s[5]] <- s[5]
    
    tmp.x <- tmp[, 4]
    s <- boxplot.stats(tmp.x)$stats
    tmp[, 4][tmp.x < s[1]] <- s[1]
    tmp[, 4][tmp.x > s[5]] <- s[5]
    out <- data.frame(cl =x, 
                      dim_1 = mean(tmp[, 3]), 
                      dim_2 = mean(tmp[, 4]), stringsAsFactors = F)
    return(out)
    
  })
  cl_coord <- do.call(rbind, cl_coord)
  
  edge.coord <- data.frame(cl1 = cl_res$S1_name,
                           cl2 = cl_res$S2_name,
                           x = cl_coord[match(cl_res$S1_name, cl_coord$cl), 2],
                           y = cl_coord[match(cl_res$S1_name, cl_coord$cl), 3],
                           xend = cl_coord[match(cl_res$S2_name, cl_coord$cl), 2],
                           yend = cl_coord[match(cl_res$S2_name, cl_coord$cl), 3])
  edge.coord <- merge(edge.coord, cl_res[, c("S1_name", "S2_name", "s")], by = c(1,2))
  edge.coord$s_norm <- (e_scale * edge.coord$s)/max(edge.coord$s)
  
  vp = unlist(lapply(1:nrow(cl_res), function(x) x <- cl_res[x, 1:2]))
  net <- graph_from_data_frame(cl_res, directed = F)
  edge.coord$edge.id = get.edge.ids(net, vp = vp) 
  dim_name_1 <- colnames(cell_emb)[3]
  dim_name_2 <- colnames(cell_emb)[4]
  
  p1 <- ggplot(cell_emb, aes( x = .data[[dim_name_1]], y = .data[[dim_name_2]], color = cluster )) +
    geom_point(size = point_size, show.legend = FALSE, alpha = point_alpha, stroke=NA) +
    geom_edge_link(data = edge.coord, aes( x = x, y = y, xend = xend, yend = yend, edge_width = s_norm), color = edge_color, alpha = edge_alpha) +
    geom_node_point(data = cl_coord, aes( x = dim_1, y = dim_2, fill = cl), shape = 21, size = node_size, show.legend = FALSE) +
    
    guides(edge_width = "none") +
    theme_test()
  
  if(!is.null(palette)) {
    p1 <- p1 + 
      scale_color_manual(values = palette, breaks = names(palette)) +
      scale_fill_manual(values = palette, breaks = names(palette))
      
      
  }
  
  if(label) {
    p1 <- p1 + geom_node_text(data = cl_coord, aes( x = dim_1, y = dim_2, label = cl),  colour = text_color, size=text_size,
                              show.legend = FALSE)
  }
  
  if(!is.null(save_file)) {
    jpeg(save_file, height = heigh, width = width, res = res, units = unit)
    print(p1)
    dev.off()
  } else {
    print(p1)
  }
  return(p1)
}


