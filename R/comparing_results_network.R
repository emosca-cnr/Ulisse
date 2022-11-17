#' Function to plot the union network of multiple PCT or CCC analyses by highlighting the shared edges and vertices.
#' @description The function is used to visually compare the results obtained by multiple PCR or CCC analyses by creating a union network and highlighting the
#' shared edges
#' @details The function takes as an input a list of data.frame resulting from multiple `pathway_cross_talk()` or `cluster_communication()` analyses
#' and produce a union network. This network is plotted by coloring the edges corresponding to how many analyses are present. The user may also decide to color
#' the vertices on the number of samples in which are present or by a user-provided variable
#' @param res_list list of resulting data.frame from `pathway_cross_talk()` or `cluster_communication()`
#' @param file_out filename of the produced plot. If not provided the plot is saved as "union_ct_network.jpeg"
#' @param colors colors to be used for the vertices. If `vertex_col_by = number` should be a vector of two colors. If `vertex_col_by` is a vector,
#' then should be a vector of colors, one for each unique variable in the `vertex_col_by` vector. If set to `NULL` then the function uses blue to red for `number`
#' and rainbow palette for vector
#' @param value value to be used for the gradient of the vertices when `vertex_col_by = number`. If `NULL` the function uses min and max
#' @param legend_vertex logical, if to print the legend of the color on the vertices
#' @param edge_width logical, if the edge width should change according to how much it is shared in the results list.
#' @param legend_col logical, if to plot the legend of the color of the edges and vertices
#' @param legend_edge_width logical, if to plot the legend of the different interval used for edge width plotting
#' @param layout layout to be used in igraph plotting. If not provided, the function calculates it
#' @param width,height,res,units graphical value of `jpeg()` function
#' @param x_vertex,y_vertex,x_edge_w,y_edge_w x and y coordinates for `legend_vertex` (when a vector is passed to `vertex_col_by` ) and `legend_edge_breaks` 
#' positioning in the plot. See `legend()` function for more details
#' @param xl,yb,xr,yt coordinates for the positioning of the gradient legend for `legend_col`. See `color.legend()` 
#' and `rect()` functions for more details
#' @param ... further argument passed to `igraph.plot` function
#' @return the function produce the plot saved with the name passed to `file_out` and also returns the layout, and the igraph object used for plotting 
#' @importFrom graphics legend
#' @importFrom grDevices rainbow dev.off jpeg
#' @importFrom plotrix color.legend
#' @importFrom circlize colorRamp2
#' @import igraph
#' @export




comparing_results_network <- function(res_list, file_out, colors =NULL, value = NULL, 
                                      edge_width =T, legend_col=T, 
                                      legend_edge_width = T, layout = NULL,
                                      width = 200, height = 200, res = 300, units = "mm", 
                                      x_vertex = -1.5, y_vertex = 0, 
                                      x_edge_w = -1.5, y_edge_w = -0.3, 
                                      xl = -1.45, yb = -0.9, xr = -1.3, yt = -0.6, ...) {
  
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
  
  Vattr.union.g <- data.frame(vertex.attributes(union.g), stringsAsFactors = F)
  
  idx <- grep("sample", colnames(Vattr.union.g))
  Vattr.union.g$nsample <- rowSums(Vattr.union.g[, idx], na.rm = T)
  V(union.g)$nsample <- Vattr.union.g$nsample
  
  #building defaults and value to be plotted---------------
  #filename
  if(is.null(file_out)) {
    file_out <- "union_ct_network.jpeg"
  }
  #vertex & edge color
  if(is.null(colors)) {
    colors <- c("blue", "red")
  }
  
  if(is.null(value)) {
    value <- c(1, length(res_list))
  }
  
  cl_fun <- colorRamp2(value, colors)
  V(union.g)$color <- cl_fun(V(union.g)$nsample)
  
  E(union.g)$color <- cl_fun(E(union.g)$nsample)
  
  
  #layout
  if(is.null(layout)) {
    lo_cl <- layout_with_fr(union.g, weights = E(union.g)$nsample)
  } else {
    lo_cl <- layout
  }
  
  #plotting----------------------------
  
  if(!edge_width) {
    jpeg(file_out, width = width, height = height, res=res, 
         units=units)
    plot(union.g, edge.col = E(union.g)$color, vertex.color=V(union.g)$color,
         layout=lo_cl, ...)
    
    if(legend_col) {
      color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                   legend = c(value[1], (value[1] + value[2])/2, value[2]), 
                   rect.col = cl_fun(seq(value[1], value[2], length.out = 10 )), 
                   gradient="y", align = "rb")
    }
    dev.off()
    
  } else {
    jpeg(file_out, width = width, height = height, res=res, 
         units=units)
    plot(union.g, edge.col = E(union.g)$color, vertex.color=V(union.g)$color,
         edge.width = (2^ E(union.g)$nsample), layout=lo_cl, ...)
    
    if(legend_edge_width) {
      legend(x= x_edge_w, y = y_edge_w, col = "black", lwd = sort(unique(2+E(union.g)$nsample)),
             legend = sort(unique(E(union.g)$nsample)), box.col = NA)
    }
    
    if(legend_col) {
      color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                   legend = c(value[1], (value[1] + value[2])/2, value[2]), 
                   rect.col = cl_fun(seq(value[1], value[2], length.out = 10 )), 
                   gradient="y", align = "rb")
    }
    
    
    dev.off()
  }
  
  return(list(network = union.g, layout = lo_cl))
  
}


