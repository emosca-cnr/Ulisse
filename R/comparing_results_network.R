#' Function to plot the union network of multiple PCT or CCC analyses by highlighting the shared edges and vertices.
#' @description The function is used to visually compare the results obtained by multiple PCR or CCC analyses by creating a union network and highlighting the
#' shared edges
#' @details The function takes as an input a list of data.frame resulting from multiple `pathway_cross_talk()` or `cluster_communication()` analyses
#' and produce a union network. This network is plotted by coloring the edges corresponding to how many analyses are present. The user may also decide to color
#' the vertices on the number of samples in which are present or by a user-provided variable
#' @param res_list list of resulting data.frame from `pathway_cross_talk()` or `cluster_communication()`
#' @param file_out filename of the produced plot. If not provided the plot is saved as "union_ct_network.jpeg"
#' @param vertex_col_by how to color the vertices. If set to `number` then the vertices are colored by the number of analyses in which is present. Otherwise 
#' can be a named vector with a value for each pathway or cell in the results. If `NULL` no color vertices will be passed to `igraph.plot` function
#' @param vertex_color colors to be sued for the vertices. If `vertex_col_by = number` should be a vector of two colors. If `vertex_col_by` is a vector,
#' then should be a vector of colors, one for each unique variable in the `vertex_col_by` vector. If set to `NULL` then the function uses blue to red for `number`
#' and rainbow palette for vector
#' @param vertex_value value to be used for the gradient of the vertices when `vertex_col_by = number`. If `NULL` the function uses min and max
#' @param legend_vertex logical, if to print the legend of the color on the vertices
#' @param edge_colors colors t be used for the edges. It should be a vector of two colors. If `NULL` the function uses blue to red
#' @param edge_value vector of the two values to be used to create the gradient to color the edges. If not provided, 
#' the function uses minimum and maximum
#' @param edge_breaks number of breaks to be used to group values in `edge_col_by` and plot the edges with different 
#' width. If set to 1, all the edges will have the same width
#' @param legend_edge_col logical, if to plot the legend of the color of the edges
#' @param legend_edge_breaks logical, if to plot the legend of the different interval used for edge width plotting
#' @param layout layout to be used in igraph plotting. If not provided, the function calculates it
#' @param width,height,res,units graphical value of `jpeg()` function
#' @param x_vertex,y_vertex,x_edge_b,y_edge_b x and y coordinates for `legend_vertex` (when a vector is passed to `vertex_col_by` ) and `legend_edge_breaks` 
#' positioning in the plot. See `legend()` function for more details
#' @param xl_vertex,yb_vertex,xr_vertex,yt_vertex,xl_edge,yb_edge,xr_edge,yt_edge coordinates for the positioning of the gradient legend for `legend_edge_col`
#' and `legend_vertex` (when `vertex_col_by = "number"` ). See `color.legend()` and `rect()` functions for more details
#' @param ... further argument passed to `igraph.plot` function
#' @return the function produce the plot saved with the name passed to `file_out` and also returns the layout, and the igraph object used for plotting 
#' @importFrom graphics legend
#' @importFrom grDevices rainbow dev.off jpeg
#' @importFrom plotrix color.legend
#' @importFrom circlize colorRamp2
#' @import igraph
#' @export




comparing_results_network <- function(res_list, file_out, vertex_col_by = NULL, vertex_colors =NULL, vertex_value = NULL, legend_vertex=F,
                                      edge_colors = NULL, edge_value = NULL, edge_breaks =1, legend_edge_col=T, 
                                      legend_edge_breaks = T, layout = NULL,
                                      width = 200, height = 200, res = 300, units = "mm", 
                                      x_vertex = -1.5, y_vertex = 0, 
                                      x_edge_b = -1.5, y_edge_b = -0.5, 
                                      xl_vertex = -1.45, yb_vertex = -0.2, xr_vertex = -1.3, yt_vertex =-0.5,
                                      xl_edge = -1.45, yb_edge = -0.9, xr_edge = -1.3, yt_edge = -0.6, ...) {
  
  res_list <- lapply(res_list, function(x) {
    x$sample <- 1
    x <- graph_from_data_frame(x, directed = F)
    V(x)$sample <- 1
    return(x)
  }
  )
  
  union.g <- igraph::union(res_list[[1]], res_list[2:length(res_list)])
  attr.union.g <- data.frame(edge.attributes(union.g), stringsAsFactors = F)
  
  idx <- grep("sample", colnames(attr.union.g))
  attr.union.g$nsample <- rowSums(attr.union.g[, idx], na.rm = T)
  
  E(union.g)$nsample <- attr.union.g$nsample
  
  
  #building defaults and value to be plotted---------------
  #filename
  if(is.null(file_out)) {
    file_out <- "union_ct_network.jpeg"
  }
  #vertex color
  if(is.null(vertex_col_by)) {
    vertex_col_by <- "number" 
  } 
  if(vertex_col_by == "number") {
    Vattr.union.g <- data.frame(vertex.attributes(union.g), stringsAsFactors = F)
    
    idx <- grep("sample", colnames(Vattr.union.g))
    Vattr.union.g$nsample <- rowSums(Vattr.union.g[, idx], na.rm = T)
    
    V(union.g)$col_by <- Vattr.union.g$nsample
    if(is.null(vertex_value)) {
      
      vertex_value <- c(min(V(union.g)$col_by), max(V(union.g)$col_by))
    }
    
    if(is.null(vertex_colors)) {
      vertex_colors <- c("blue", "red")
    }
    if(vertex_value[1] == vertex_value[2]) {
      V(union.g)$color <- "red"
    } else {
      vertex_cl <- colorRamp2(vertex_value, vertex_colors)
      V(union.g)$color <- vertex_cl(V(union.g)$col_by)
    }
    
  } else {
    V(union.g)$col_by <- vertex_col_by[V(union.g)$name]
    if(!is.null(vertex_color)) {
      if(length(V(union.g)$name) > length(vertex_col_by)) {
        print.warnings("Not enought colors provided for the vertices. Rainbow used instead")
        vertex_color <- rainbow(length(V(union.g)$name))
        names(vertex_color) <- unique(V(union.g)$name)
      } 
      names(vertex_color) <- unique(V(union.g)$name)
    } else {
      vertex_color <- rainbow(length(V(union.g)$name))
      names(vertex_color) <- unique(V(union.g)$name)
    }
    V(union.g)$color <- vertex_color[V(union.g)$name]
  }
  #edges color
  if(is.null(edge_colors)) {
    edge_colors <- c("blue", "red")
  }
  
  if(is.null(edge_value)) {
    edge_value <- c(min(E(union.g)$nsample), max(E(union.g)$nsample))
  }
  
  edge_cl <- colorRamp2(edge_value, edge_colors)
  E(union.g)$color <- edge_cl(E(union.g)$nsample)
  
  if(edge_breaks==1) {
    edge_breaks <- NULL
    } else {
      edge_breaks <- cut(E(union.g)$nsample, breaks = edge_breaks)
    }
  
  #layout
  if(is.null(layout)) {
    lo_cl <- layout_with_fr(union.g, weights = E(union.g)$nsample)
  } else {
    lo_cl <- layout
  }
  
  #plotting----------------------------
  if(is.null(vertex_col_by)){
    if(is.null(edge_breaks)) {
      jpeg(file_out, width = width, height = height, res=res, 
           units=units)
      plot(union.g, edge.col = E(union.g)$color, layout=lo_cl, ...)
      if(legend_edge_col) {
        color.legend(xl = xl_edge, yb = yb_edge, xr = xr_edge, yt = yt_edge, 
                     legend = c(edge_value[1], (edge_value[1] + edge_value[2])/2, edge_value[2]), 
                     rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                     gradient="y", align = "rb")
      }
      
      dev.off()
    } else {
      jpeg(file_out, width = width, height = height, res=res, 
           units=units)
      plot(union.g, edge.col = E(union.g)$color, 
           edge.width = (2+ as.numeric(edge_breaks)), layout=lo_cl, ...)
      if(legend_edge_breaks) {
        legend(x= x_edge_b, y = y_edge_, col = "black", lwd = sort(unique(2+as.numeric(edge_breaks))),
               legend = levels(edge_breaks), box.col = NA)
      }
      if(legend_edge_col) {
        color.legend(xl = xl_edge, yb = yb_edge, xr = xr_edge, yt = yt_edge, 
                     legend = c(edge_value[1], (edge_value[1] + edge_value[2])/2, edge_value[2]), 
                     rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                     gradient="y", align = "rb")
      }
      dev.off()
      
    }
  } else {
    if(vertex_col_by == "number") {
      if(is.null(edge_breaks)) {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(union.g, edge.col = E(union.g)$color, vertex.color=V(union.g)$color,
             layout=lo_cl, ...)
        if(legend_vertex) {
          color.legend(xl = xl_vertex, yb = yb_vertex, xr = xr_vertex, yt = yt_vertex, 
                       legend = c(vertex_value[1], (vertex_value[1] + vertex_value[2])/2, vertex_value[2]), 
                       rect.col = edge_cl(seq(vertex_value[1], vertex_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        
        if(legend_edge_col) {
          color.legend(xl = xl_edge, yb = yb_edge, xr = xr_edge, yt = yt_edge, 
                       legend = c(edge_value[1], (edge_value[1] + edge_value[2])/2, edge_value[2]), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        dev.off()
        
      } else {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(union.g, edge.col = E(union.g)$color, vertex.color=V(union.g)$color,
             edge.width = (2+ as.numeric(edge_breaks)), layout=lo_cl, ...)
        
        if(legend_edge_breaks) {
          legend(x= x_edge_b, y = x_edge_b, col = "black", lwd = sort(unique(2+as.numeric(edge_breaks))),
                 legend = levels(edge_breaks), box.col = NA)
        }
        
        if(legend_vertex) {
          color.legend(xl = xl_vertex, yb = yb_vertex, xr = xr_vertex, yt = yt_vertex, 
                       legend = c(vertex_value[1], (vertex_value[1] + vertex_value[2])/2, vertex_value[2]), 
                       rect.col = edge_cl(seq(vertex_value[1], vertex_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        if(legend_edge_col) {
          color.legend(xl = xl_edge, yb = yb_edge, xr = xr_edge, yt = yt_edge, 
                       legend = c(edge_value[1], (edge_value[1] + edge_value[2])/2, edge_value[2]), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        
        
        dev.off()
      }
    } else {
      if(is.null(edge_breaks)) {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(union.g, edge.col = E(union.g)$color, vertex.color=V(union.g)$color,
             layout=lo_cl, ...)
        if(legend_vertex) {
          legend(x= x_vertex, y = x_vertex, fill = vertex_colors, legend = names(vertex_colors))
        }
        
        if(legend_edge_col) {
          color.legend(xl = xl_edge, yb = yb_edge, xr = xr_edge, yt = yt_edge, 
                       legend = c(edge_value[1], (edge_value[1] + edge_value[2])/2, edge_value[2]), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        dev.off()
      } else {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(union.g, edge.col = E(union.g)$color, vertex.color=V(union.g)$color,
             edge.width = (2+ as.numeric(edge_breaks)), layout=lo_cl, ...)
        
        if(legend_edge_breaks) {
          legend(x= x_edge_b, y = x_edge_b, col = "black", lwd = sort(unique(2+as.numeric(edge_breaks))),
                 legend = levels(edge_breaks), box.col = NA)
        }
        if(legend_vertex) {
          legend(x= x_vertex, y = x_vertex, fill = vertex_colors, legend = names(vertex_colors))
        }
        
        if(legend_edge_col) {
          color.legend(xl = xl_edge, yb = yb_edge, xr = xr_edge, yt = yt_edge, 
                       legend = c(edge_value[1], (edge_value[1] + edge_value[2])/2, edge_value[2]), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        
        
        dev.off()
      }
    }
    
    
  }
  
  return(list(network = union.g, layout = lo_cl))
  
}


