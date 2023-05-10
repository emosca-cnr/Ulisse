#' Function to plot the heatmap of cross-talk results
#' @description The function elaborates the output of `gs_cross-talk()` 
#' function to represent the obtained cross-talks as an heatmap
#' @details The functions uses `gs_cross-talk()` output to build a cross-talk network. Then, 
#' the adjacency matrix of the network is plotted with or without annotations on the rows and columns. The user may 
#' also decide which value to plot in the heatmap
#' @param ct results obtained with `gs_cross-talk()` functions
#' @param color_by name of the column to be plotted in the heatmap
#' @param color vector of two colors used to build the gradient to color the heatmap, that are the colors of the two
#' number of `color_level`. If NULL, the function uses `lightyellow` and `red3`
#' @param no_ct_color color that have to be used for the gene set pairs that shows no cross-talk (in the heatmap will have score = 0). If `NULL` the function
#' uses `whitesmoke`
#' @param color_level vector of the two values to be used to create the gradient to color the heatmap. If not provided, 
#' the function uses minimum and maximum
#' @param filtering logical, if the function have to plot all the edges or only the significant ones. If `TRUE`, the function uses `p_val, FDR` and `ct_val` to identify the 
#'  significant ones, that are plotted with a solid line, and the not significant with a dashed line. If `FALSE`, all edges are plotted as a solid line. 
#' @param p_val,FDR,ct_val filtering values. If one of these is set to `NULL` the function ignores it
#' @param legend_side where to place the legend of the heatmap and of the annotations (if provided). See `Heatmap` 
#' for further details
#' @param community logical or an object resulting from igraph community calculation. If `TRUE`, the function calculates
#' the communities by using fastgreedy algorithm, if `FALSE` the function does not consider communities. Otherwise, if a
#' community object is provided, the function uses it and colors the vertices accordingly. The communities will be
#' plotted as a row and column annotation of the heatmap
#' @param pal_community vector of colors to be used to color the community annotation. If `NULL` the function will use 
#' `rainbow()` palette
#' @param label_size size of the gene-set names that will be printed on the diagonal of the heatmaps
#' @param row_annotation,column_annotation data.frame with columns corresponding to the rows or columns annotations. The rownames must 
#' be named after all the cells/gene-set in the `ct` tables. The values in the columns will be treated as a 
#' discrete variable
## #' @param row_annotation_name,column_annotation_name name to be used as a title for the legend of the row/column annotation
#' @param pal_row_annotation,pal_column_annotation list of vectors, one for each column in `row_annotation,column_annotation` df. Each 
#'  vector should be a color vector named after each cells/gene-set in the `ct` table. If `NULL`, the function uses `pals::alphabet2()` palette
#' @param row_name_side,col_name_side if the name of the row/column annotation should be put at the "bottom" or "top" ("left" or "right") 
#' of the annotation
#' @param file_out name used to save the jpeg file. If `NULL` the complexHeatmap object is returned
#' @param file_width,file_height,res,units graphical value of `jpeg()` function
#' @param ... further arguments passed down to `ComplexHeatmap::Heatmap()`
#' @return the function produce the plot saved with the name passed to `file_out` and also returns the adjacency matrix 
#' and the communities (if calculated)
#' @importFrom grDevices rainbow dev.off jpeg
#' @importFrom circlize colorRamp2
#' @import igraph
#' @import ComplexHeatmap
#' @import grid
#' @export


ct_heatmap <- function(ct, color_by = "ct_score", color = NULL, color_level = NULL, no_ct_color = NULL,
                       filtering = FALSE, p_val, FDR, ct_val,
                       legend_side = "left", community = F, pal_community = NULL, label_size = 5,
                       row_annotation = NULL, pal_row_annotation = NULL, row_name_side = "bottom",
                       column_annotation = NULL, pal_column_annotation = NULL, col_name_side = "left",
                       file_out = NULL, file_width = 200, file_height = 200, res = 300, units = "mm", ...) {
  
  ct.g <- graph_from_data_frame(ct, directed = F)
  adj.ct <- as_adjacency_matrix(ct.g, sparse = F, attr = color_by)
  ct_names <- as.character(rownames(adj.ct))
  if(filtering) {
    if(is.null(p_val)) {
      p_val <- max(ct$p_value_link, na.rm = T)
    }
    if(is.null(FDR)) {
      FDR <- max(ct$FDR_link, na.rm = T)
    }
    if(is.null(ct_val)) {
      ct_val <- min(ct$ct_score, na.rm = T)
    }
    
    sig.list <- list(adj = sign(as_adjacency_matrix(ct.g, attr = "ct_score")),
                     p.val = as_adjacency_matrix(ct.g, attr = "p_value_link" ),
                     eFDR = as_adjacency_matrix(ct.g, attr = "FDR_link"))
    
    sig.ct <- matrix(data = 0, nrow = length(ct_names), ncol = ncol(adj.ct), dimnames = list(ct_names, ct_names))
    for(i in 1:length(ct_names)) {
      sig.ct[i,][sig.list$p.val[i,] <= p_val & sig.list$eFDR[i,] <= FDR & sig.list$adj[i,] != 0 & adj.ct[i,] >= ct_val] <- 1
      
    }
  } else {
    sig.ct <- matrix(data = 0, nrow = length(ct_names), ncol = ncol(adj.ct), dimnames = list(ct_names, ct_names))
  }
  
  #setting options--------------------
  #communities
  if(is.logical(community)) {
    if(community) {
      community <- fastgreedy.community(ct.g, weights = edge_attr(ct.g)[["ct_score"]])
    } else {
      community <- NULL
    }
  }
  if(!is.null(community)) {
    community <- membership(community)
    comm_ct_mem <- as.factor(community)
    names(comm_ct_mem) <- V(ct.g)$name
  }
  
  if(is.null(pal_community) & !is.null(community)) {
    pal_community <- rainbow(n = max(as.numeric(levels(comm_ct_mem))))
    names(pal_community) <- 1:max(as.vector(as.numeric(community)))
  } else if(!is.null(pal_community) & !is.null(community)) {
    n <- max(as.numeric(community))
    if(length(pal_community) < n) {
      print.warnings("Not enought color provided in the community palette")
    } else {
      pal_community <- pal_community[1:n]
      names(pal_community) <- 1:max(as.vector(as.numeric(community)))
    }
  }
  
  #row annotation----------------------
  if(!is.null(row_annotation)) {
    df.row <- row_annotation[match(rownames(adj.ct), rownames(row_annotation)),, drop = F]
    if(is.null(pal_row_annotation)) {
      pal_row <- lapply(1:ncol(row_annotation), function(x) {
        n <- unique(row_annotation[, x])
        pal <- pals::alphabet2(length(n))
        names(pal) <- n
        return(pal)
      })
      names(pal_row) <- colnames(row_annotation)
    }
    
    if(!is.null(community)) {
      df.row$community <- comm_ct_mem[df.row$name]
      pal_row <- c(pal_row, community = list(pal_community))
    }
  } else {
    if(!is.null(community)) {
      comm_ct_mem <- sort(comm_ct_mem)
      adj.ct <- adj.ct[names(comm_ct_mem), names(comm_ct_mem), drop = F]
      df.row <- data.frame(community = comm_ct_mem, stringsAsFactors = F)
      pal_row <- list(community = pal_community)
    }
    
  }
  
  
  #column annotation---------------
  if(!is.null(column_annotation)) {
    df.col <- column_annotation[match(rownames(adj.ct), rownames(column_annotation)),, drop = F]
    if(is.null(pal_column_annotation)) {
      pal_column <- lapply(1:ncol(row_annotation), function(x) {
        n <- unique(row_annotation[, x])
        pal <- pals::alphabet2(length(n))
        names(pal) <- n
        return(pal)
      })
      names(pal_column) <- colnames(column_annotation)
    }
    
    if(!is.null(community)) {
      df.col$community <- comm_ct_mem[df.col$name]
      pal_column <- c(pal_column, community = list(pal_community))
    }
  } else {
    if(!is.null(community)) {
      comm_ct_mem <- sort(comm_ct_mem)
      adj.ct <- adj.ct[names(comm_ct_mem), names(comm_ct_mem), drop = F]
      df.col <- data.frame(community = comm_ct_mem)
      pal_column<- list(community = pal_community)
    }
    
  }
  
  #file name 
  # if(is.null(file_out)) {
  #   file_out <- "CT_heatmap.jpeg"
  # }
  #color function
  if(is.null(color)) {
    color <- c("blue", "red")
  }
  
  if(!is.null(color_by)) {
    edge_att <- edge_attr(ct.g, name = color_by)
    if(is.null(color_level)) {
      
      color_level <- c(min(edge_att), max(edge_att))
    }
  }
  
  if(is.null(no_ct_color)) {
    no_ct_color <- "whitesmoke"
  }
  col_fun <- colorRamp2(c(0, color_level), c(no_ct_color, color))
  
  #plotting--------------------
  if(!is.null(file_out)) {
    jpeg(file_out, width = file_width, height = file_height, res=res, 
         units=units)
  }
  
  
  if(is.null(column_annotation)) {
    if(is.null(row_annotation)) {
      
      h <- Heatmap(adj.ct, col = col_fun, rect_gp = gpar(type = "none"),
                   heatmap_legend_param = list(title=color_by), cell_fun = function(j, i, x, y, width, height, fill) {
                     
                     if(i > j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = "white", fill = col_fun(adj.ct[i, j]), lwd = 0.2))
                       if(sig.ct[i, j] == 1) {
                         grid.text("*", x, y, gp = gpar(col="blue"))
                       }
                       
                     } else if(i < j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = NA,  fill = NA))
                     } else {
                       txt <- rownames(adj.ct)[i]
                       grid.text(txt, x, y, gp = gpar(fontsize = label_size), just = "left")
                     }
                   }, cluster_rows = FALSE, cluster_columns = FALSE,
                   show_row_names = FALSE, show_column_names = FALSE, ...)
      h <- draw(h, heatmap_legend_side = legend_side, merge_legend =T)
      
    } else {
      row_ha <- HeatmapAnnotation(df = df.row, col = pal_row, annotation_name_side = row_name_side, which = "row")
      
      h <- Heatmap(adj.ct, col = col_fun,
                   heatmap_legend_param = list(title=color_by),
                   left_annotation = row_ha, rect_gp = gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
                     
                     if(i > j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = "white", fill = col_fun(adj.ct[i, j]), lwd = 0.2))
                       if(sig.ct[i, j] == 1) {
                         grid.text("*", x, y, gp = gpar(col="blue"))
                       }
                       
                     } else if(i < j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = NA,  fill = NA))
                     } else {
                       txt <- rownames(adj.ct)[i]
                       grid.text(txt, x, y, gp = gpar(fontsize = label_size), just = "left")
                     }
                   }, cluster_rows = FALSE, cluster_columns = FALSE,
                   show_row_names = FALSE, show_column_names = FALSE, ...)
      h <- draw(h, heatmap_legend_side = legend_side, merge_legend =T)
    }
  } else {
    if(is.null(row_annotation)) {
      
      col_ha <- columnAnnotation(df = df.col, col = pal_column, annotation_name_side = col_name_side)
      h <- Heatmap(adj.ct, col = col_fun,
                   heatmap_legend_param = list(title=color_by),
                   bottom_annotation = col_ha,
                   rect_gp = gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
                     
                     if(i > j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = "white", fill = col_fun(adj.ct[i, j]), lwd = 0.2))
                       if(sig.ct[i, j] == 1) {
                         grid.text("*", x, y, gp = gpar(col="blue"))
                       }
                       
                     } else if(i < j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = NA,  fill = NA))
                     } else {
                       txt <- rownames(adj.ct)[i]
                       grid.text(txt, x, y, gp = gpar(fontsize = label_size), just = "left")
                     }
                   }, cluster_rows = FALSE, cluster_columns = FALSE,
                   show_row_names = FALSE, show_column_names = FALSE, ...)
      h <- draw(h, heatmap_legend_side = legend_side, merge_legend =T)
      
    } else {
      
      row_ha <- HeatmapAnnotation(df = df.row, col = pal_row, annotation_name_side = row_name_side, which = "row")
      col_ha <- HeatmapAnnotation(df = df.col, col = pal_column, annotation_name_side = col_name_side)
      h <- Heatmap(adj.ct, col = col_fun,
                   heatmap_legend_param = list(title=color_by),
                   left_annotation = row_ha, 
                   bottom_annotation = col_ha,
                   rect_gp = gpar(type = "none"), cell_fun = function(j, i, x, y, width, height, fill) {
                     
                     if(i > j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = "white", fill = col_fun(adj.ct[i, j]), lwd = 0.2))
                       if(sig.ct[i, j] == 1) {
                         grid.text("*", x, y, gp = gpar(col="blue"))
                       }
                       
                     } else if(i < j) {
                       grid.rect(x = x, y = y, width = width, height = height,
                                 gp = gpar(col = NA,  fill = NA))
                     } else {
                       txt <- rownames(adj.ct)[i]
                       grid.text(txt, x, y, gp = gpar(fontsize = label_size), just = "left")
                     }
                   }, cluster_rows = FALSE, cluster_columns = FALSE,
                   show_row_names = FALSE, show_column_names = FALSE)#, ...)
      h <- draw(h, heatmap_legend_side = legend_side, merge_legend =T)
      
    }
  }
  
  print(h)
  if(!is.null(file_out)) {
    dev.off()
  }
  
  if(is.null(file_out)) {
    return(h)
  }
}
