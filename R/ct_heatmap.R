#' Function to plot the heatmap of CT network
#' @description The function elaborate the output of `pathway_cross_talk()` or `cluster_communication()` 
#' function to represent the cross-talk network obtained as an heatmap
#' @details The functions uses `pathway_cross_talk()` or `cluster_communication()` output to build a CT network. Then, 
#' the adjacency matrix of the network is plotted with or without annotations on the rows and columns. The user may 
#' also decide which value to plot in the heatmap
#' @param ct output of `pathway_cross_talk()` or `cluster_communication()` function
#' @param color_by name of the column to be plotted in the heatmap
#' @param color vector of two colors used to build the gradient to color the heatmap. If NULL, the function uses
#' "blue" to "red" 
#' @param color_level vector of the two values to be used to create the gradient to color the heatmap. If not provided, 
#' the function uses minimum and maximum
#' @param legend_side where to place the legend of the heatmap and of the annotations (if provided). See `Heatmap` 
#' for further details
#' @param community logical or an object resulting from igraph community calculation. If TRUE, the function calculates
#' the communities by using fastgreedy algorithm, if FALSE the function does not consider communities. Otherwise, if a
#' community object is provided, the function use it and color the vertices accordingly. The communities will be
#' plotted as a row and column annotation of the heatmap
#' @param pal_community vector of colors to be used to color the community annotation. If NULL the function will use 
#' `rainbow()` palette
#' @param row_annotation,column_annotation named vector to be used to annotate the rows or columns. The vector must 
#' be named after all the cells/pathways in the `ct` tables. The values in the vector will be treated as a 
#' discrete variable
#' @param row_annotation_name,column_annotation_name name to be used as a title for the legend of the row/column annotation
#' @param pal_row_annotation,pal_column_annotation vector of colors to be used to color the annotation. If NULL, the 
#' function uses `rainbow()` palette
#' @param row_name_side,col_name_side if the name of the row(column) annotation should be put at the "bottom" or "top" ("left" or "right") 
#' of the annotation
#' @param file_out name used to save the jpeg file. If not provided, the plot will be saved as "CT_heatmap.jpeg"
#' @param width,height,res,units graphical value of `jpeg()` function
#' @param ... further arguments passed down to ComplexHeatmap::Heatmap
#' @return the function produce the plot saved with the name passed to `file_out` and also returns the adjacency matrix 
#' and the communities (if calculated)
#' @importFrom grDevices rainbow dev.off jpeg
#' @importFrom circlize colorRamp2
#' @import igraph
#' @import ComplexHeatmap
#' @export


ct_heatmap <- function(ct, color_by = "pct", color = NULL, color_level = NULL, legend_side = "left",
                       community = F, pal_community = NULL, 
                       row_annotation = NULL, row_annotation_name = "row_anno", pal_row_annotation = NULL, row_name_side = "bottom",
                       column_annotation = NULL, column_annotation_name = "col_anno", pal_column_annotation = NULL, col_name_side = "left",
                       file_out = NULL, width = 200, height = 200, res = 300, units = "mm", ...) {
  
  ct.g <- graph_from_data_frame(ct, directed = F)
  adj.ct <- as_adjacency_matrix(ct.g, sparse = F, attr = color_by)
  
  ct_column <- grep("pct|ccc_score", colnames(ct))
  ct_column <- colnames(ct)[ct_column]
  #setting options--------------------
  #communities
  if(is.logical(community)) {
    if(community) {
      community <- fastgreedy.community(ct.g, weights = edge_attr(ct.g)[[ct_column]])
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
  
  #row annotation
  if(!is.null(row_annotation)) {
    df.row <- data.frame(names(row_annotation), row_annotation, stringsAsFactors = F)
    colnames(df.row) <- c("name", row_annotation_name)
    if(is.null(pal_row_annotation)) {
      pal_row_annotation <- rainbow(length(unique(row_annotation)))
    }
    names(pal_row_annotation) <- unique(names(row_annotation))
    pal_row <- list(pal_row_annotation)
    names(pal_row) <- row_annotation_name
    if(!is.null(community)) {
      df.row$community <- comm_ct_mem[df.row$name]
      pal_row <- c(pal_row, list(pal_community))
      names(pal_row) <- c(row_annotation_name, "community")
    }
  } else {
    if(!is.null(community)) {
      comm_ct_mem <- sort(comm_ct_mem)
      adj.ct <- adj.ct[names(comm_ct_mem), names(comm_ct_mem)]
      df.row <- data.frame(name = names(comm_ct_mem), community = comm_ct_mem)
      pal_row <- list(community = pal_community)
      
    }
    
  }
  
  
  #column annotation
  if(!is.null(column_annotation)) {
    df.col <- data.frame(names(column_annotation), column_annotation, stringsAsFactors = F)
    colnames(df.row) <- c("name", column_annotation_name)
    if(is.null(pal_column_annotation)) {
      pal_column_annotation <- rainbow(length(unique(column_annotation)))
    }
    names(pal_column_annotation) <- unique(names(column_annotation))
    pal_column <- list(pal_column_annotation)
    names(pal_column) <- column_annotation_name
    if(!is.null(community)) {
      df.col$community <- comm_ct_mem[df.col$name]
      pal_column <- c(pal_column, list(pal_community))
      names(pal_column) <- c(column_annotation_name, "community")
    }
  } else {
    if(!is.null(community)) {
      comm_ct_mem <- sort(comm_ct_mem)
      adj.ct <- adj.ct[names(comm_ct_mem), names(comm_ct_mem)]
      df.col <- data.frame(name = names(comm_ct_mem), community = comm_ct_mem)
      pal_column<- list(community = pal_community)
      
    }
    
  }
  
  #file name 
  if(is.null(file_out)) {
    file_out <- "CT_heatmap.jpeg"
  }
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
  
  col_fun <- colorRamp2(color_level, color)
  
  #plotting--------------------
  
  if(is.null(column_annotation)) {
    if(is.null(row_annotation)) {
      jpeg(file_out, width = width, height = height, res=res, 
           units=units)
      h <- Heatmap(adj.ct, col = col_fun,
                   heatmap_legend_param = list(title=color_by),
                   ...)
      draw(h, heatmap_legend_side = legend_side, merge_legend =T)
      dev.off()
    } else {
      jpeg(file_out, width = width, height = height, res=res, 
           units=units)
      row_ha <- rowAnnotation(df.row, col = pal_row, annotation_name_side = row_name_side)
      h <- Heatmap(adj.ct, col = col_fun,
                   heatmap_legend_param = list(title=color_by),
                   right_annotation = row_ha, 
                   ...)
      draw(h, heatmap_legend_side = legend_side, merge_legend =T)
      dev.off()
    }
  } else {
    if(is.null(row_annotation)) {
      jpeg(file_out, width = width, height = height, res=res, 
           units=units)
      col_ha <- columnAnnotation(df.col, col = pal_column, annotation_name_side = col_name_side)
      h <- Heatmap(adj.ct, col = col_fun,
                   heatmap_legend_param = list(title=color_by),
                   bottom_annotation = col_ha,
                   ...)
      draw(h, heatmap_legend_side = legend_side, merge_legend =T)
      dev.off()
    } else {
      jpeg(file_out, width = width, height = height, res=res, 
           units=units)
      row_ha <- rowAnnotation(df.row, col = pal_row, annotation_name_side = row_name_side)
      col_ha <- columnAnnotation(df.col, col = pal_column, annotation_name_side = col_name_side)
      h <- Heatmap(adj.ct, col = col_fun,
                   heatmap_legend_param = list(title=color_by),
                   right_annotation = row_ha, 
                   bottom_annotation = col_ha,
                   ...)
      draw(h, heatmap_legend_side = legend_side, merge_legend =T)
      dev.off()
    }
  }
  
}
#   if(is.null(community)) {
#     if(is.null(column_annotation)) {
#       if(is.null(row_annotation)) {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      ...)
#         draw(h, heatmap_legend_side = legend_side)
#         dev.off()
#       } else {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         row_ha <- rowAnnotation(df.row, col = pal_row)
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      right_annotation = row_ha, 
#                      ...)
#         draw(h, heatmap_legend_side = legend_side)
#         dev.off()
#       }
#     } else {
#       if(is.null(row_annotation)) {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         col_ha <- columnAnnotation(df.col, col = pal_column)
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      bottom_annotation = col_ha,
#                      ...)
#         draw(h, heatmap_legend_side = legend_side)
#         dev.off()
#       } else {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         row_ha <- rowAnnotation(df.row, col = pal_row)
#         col_ha <- columnAnnotation(df.col, col = pal_column)
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      right_annotation = row_ha, 
#                      bottom_annotation = col_ha,
#                      ...)
#         draw(h, heatmap_legend_side = legend_side)
#         dev.off()
#       }
#     }
#   } else {
#     if(is.null(column_annotation)) {
#       if(is.null(row_annotation)) {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         comm_ct_mem <- sort(comm_ct_mem)
#         adj.ct <- adj.ct[names(comm_ct_mem), names(comm_ct_mem)]
#         row_ha <- rowAnnotation(community = comm_ct_mem, col = list(community = pal_community), annotation_name_side = "bottom")
#         col_ha <- columnAnnotation(community = comm_ct_mem, col = list(community = pal_community), annotation_name_side = "left")
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      right_annotation = row_ha, 
#                      bottom_annotation = col_ha,
#                      ...)
#         draw(h, heatmap_legend_side = legend_side, annotation_legend_side = legend_side, merge_legend =T)
#         dev.off()
#       } else {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         row_ha <- rowAnnotation(df.row, col = pal_row)
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      right_annotation = row_ha, 
#                      ...)
#         draw(h, heatmap_legend_side = legend_side)
#         dev.off()
#        
#       }
#     } else {
#       if(is.null(row_annotation)) {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         col_ha <- columnAnnotation(df.col, col = pal_column)
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      bottom_annotation = col_ha,
#                      ...)
#         draw(h, heatmap_legend_side = legend_side)
#         dev.off()
#       } else {
#         jpeg(file_out, width = width, height = height, res=res, 
#              units=units)
#         row_ha <- rowAnnotation(df.row, col = pal_row)
#         col_ha <- columnAnnotation(df.col, col = pal_column)
#         h <- Heatmap(adj.ct, col = col_fun,
#                      heatmap_legend_param = list(title=color_by),
#                      right_annotation = row_ha, 
#                      bottom_annotation = col_ha,
#                      ...)
#         draw(h, heatmap_legend_side = legend_side)
#         dev.off()
#       }
#     }
#   }
#   if(!is.null(community)) {
#     return(list(matrix = adj.ct, community = community))
#   } else {
#     return(adj.ct)
#   }
#    
# }
# 
# 




