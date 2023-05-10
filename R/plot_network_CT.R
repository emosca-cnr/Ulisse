#' Plot cross-talk results as a network
#' @description The function elaborates the results of `gs_cross_talk()` functions to obtain a graphical 
#'  representation of the cross-talk network
#' @details The functions uses `gs_cross_talk()` output to build a cross-talk network. Then, the network is plotted with
#'  the edges colored by a value of interest. The user may decide if to plot all the edges, color the vertices by 
#'  their community, scale the dimension of the edges over a variable of interest.
#' @param ct output of `gs_cross_talk()` function
#' @param filtering logical, if the function have to plot all the edges or only the significant ones. If `TRUE`, the function uses `p_val, FDR` and `ct_val` to identify the 
#'  significant ones, that are plotted with a solid line, and the not significant with a dashed line. If `FALSE`, all edges are plotted as a solid line. 
#' @param p_val,FDR,ct_val filtering values. If one of these is set to `NULL` the function ignores it
#' @param community logical or an object resulting from igraph community calculation. If `TRUE`, the function calculates
#'  the communities by using fastgreedy algorithm, if `FALSE` the function does not consider communities. Otherwise, if a
#'  community object is provided, the function uses it and color the vertices accordingly
#' @param pal_community color palette to be used to color the vertices according to the communities. If not provided,
#'  the function uses `rainbow()` palette. If `community = FALSE` the function ignores it
#' @param vertex a list of two vector, where the first is the name of the attribute and the second is a named vector with the attribute per vertices 
#'  to be used for color (only discrete value), or `NULL` if the vertices should be not colored (thus "grey65" is passed as default color)
#' @param vertex_pal color palette used to color vertices. Should be a named vector with the names corresponding to the unique elements in `vertices`. If `NULL`
#'  the functions will use pals::alphabet2 as default
#' @param vertex_label logical, if vertex names should be plotted or not, or a named vector, with the label that should be associated to each vertex named by them
#' @param vertex_size size of the vertices in the plot
#' @param voronoi_radius,voronoi_alpha: parameters passed to `geom_node_voronoi()`. If both `vertex` and `community` arguments are enabled, the the communities are
#'  represented as cells of voronoi tessellation (see `ggraph::geom_node_voronoi()`). In this case, these argument are used to control 
#'  the radius of the cells and the transparency, respectively.
#' @param edge_col_by name of the column to be used to color the edges. Default = `ct_score`
#' @param edge_pal color to be used to create the gradient to color the edges. If not provided, the function uses "blue" 
#'  to "red"
#' @param edge_width logical, if the edge widths should be proportional to `edge_col_by` values
#' @param edge_adj_col value used to adjust color transparency of the edges
#' @param file_out name of the jpeg file produced. if `NULL` the functions returns also the plot object
#' @param width,height,res,units graphical value of `jpeg()` function
#' @param ... further graphical parameters to be passed to `ggraph()` function
#' @return If `file_out` is null the function returns the plot and the igraph object used for plotting (which may contain the communities under `comm_id` attribute). 
#'  Otherwise, only the igraph object is returned and the plot is saved through to `file_out` 
#' @importFrom grDevices rainbow dev.off jpeg adjustcolor
#' @importFrom stats setNames
#' @import igraph
#' @import ggraph
#' @export


plot_network_CT <- function(ct, filtering = FALSE, p_val, FDR, ct_val, 
                            community, pal_community =NULL, voronoi_radius = 0.8, voronoi_alpha = 0.3,
                            vertex = "name", vertex_pal = NULL, vertex_size = 5, vertex_label = TRUE,
                            edge_col_by = "ct_score", edge_pal=NULL, 
                            edge_width = T, edge_adj_col = 0.7,
                            file_out =NULL,  width = 200, height = 200, res = 300, units = "mm", ...) {
  
  #building defaults and values to be plotted--------------
  
  #filtering
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
    ct$filt <- "not"
    ct$filt[ct$p_value_link <= p_val & ct$FDR_link <= FDR & ct$ct_score >= ct_val] <- "yes"
  }
  #network--------------------------
  ct.g <- graph_from_data_frame(ct, directed = F)
  #community-------------------------
  if(is.logical(community)) {
    if(!community) {
      community <- NULL
    } else {
      community <- fastgreedy.community(ct.g, weights = edge_attr(ct.g)[["ct_score"]])
    }
  }
  
  if(!is.null(community)) {
    community <- igraph::membership(community)
    comm_ct_mem <- as.character(community)
    names(comm_ct_mem) <- V(ct.g)$name
    V(ct.g)$comm_id <- comm_ct_mem
  }
  
  if(is.null(pal_community) & !is.null(community)) {
    pal_community <- rainbow(n = max(as.numeric(unique(comm_ct_mem))))
    names(pal_community) <- sort(as.vector(as.numeric(unique(comm_ct_mem))))
    } else if(!is.null(pal_community) & !is.null(community)) {
    n <- max(as.numeric(community))
    if(length(pal_community) < n) {
      print.warnings("Not enought color provided in the community palette")
    } else {
      pal_community <- pal_community[1:n]
      names(pal_community) <- sort(as.vector(as.numeric(unique(comm_ct_mem))))
    }
    }
  #vertices---------------------
  if(!is.null(vertex)) {
  if(length(vertex) == 2) {
    vertex_attr(ct.g, name = vertex[[1]]) <- vertex[[2]][V(ct.g)$name]
    vertex <- vertex[[1]]
  }
  }
  if(!is.null(vertex) & is.null(vertex_pal)) {
    v <- unique(vertex_attr(ct.g, name = vertex[[1]]))
    vertex_pal <- setNames(pals::alphabet2(length(v)), v)
  }
  
  if(is.logical(vertex_label)) {
    if(vertex_label) {
      V(ct.g)$label <- V(ct.g)$name
    }
  } else {
    V(ct.g)$label <- vertex_label[V(ct.g)$name]
    vertex_label <- TRUE
  }
  
  #edges--------------------------
  if(is.null(edge_pal)) {
    edge_pal <- c("blue", "red")
  }
  
  #plotting-----------------------------
  
  p <- ggraph(ct.g, ...) +
    theme_graph() 
  if(edge_width & filtering) {
    p <- p + geom_edge_link(aes_string(edge_colour = edge_col_by, edge_width = edge_col_by, linetype = "filt"), alpha = edge_adj_col) +
      scale_edge_color_gradient(low = edge_pal[1], high = edge_pal[2]) +
      scale_edge_linetype_manual(values = c("dashed", "solid"), breaks = c("not", "yes"), labels = c("Not significant", "Significant"), name = "Filtering")
      
  } else if(edge_width & !filtering) {
    p <- p + geom_edge_link(aes_string(edge_colour = edge_col_by, edge_width = edge_col_by), alpha = edge_adj_col) +
      scale_edge_color_gradient(low = edge_pal[1], high = edge_pal[2])
    
  } else if(!edge_width & filtering) {
    p <- p + geom_edge_link(aes_string(edge_colour = edge_col_by, linetype = "filt"), alpha = edge_adj_col) +
      scale_edge_color_gradient(low = edge_pal[1], high = edge_pal[2]) +
      scale_edge_linetype_manual(values = c("dashed", "solid"), breaks = c("not", "yes"), labels = c("Not significant", "Significant"), name = "Filtering")
    
  } else if (!edge_width & !filtering) {
    p <- p + geom_edge_link(aes_string(edge_colour = edge_col_by), alpha = edge_adj_col) +
      scale_edge_color_gradient(low = edge_pal[1], high = edge_pal[2]) 
  }
  
  if(!is.null(community) & !is.null(vertex)) {
    p <- p + geom_node_voronoi(aes_string(fill = "comm_id"), max.radius = voronoi_radius, colour = 'white', alpha = voronoi_alpha) +
      scale_fill_manual(limits = names(pal_community), values = pal_community) +
      geom_node_point(aes_string(color = vertex), size = vertex_size) +
      scale_color_manual(limits = names(vertex_pal), values = vertex_pal)
      
  } else if(!is.null(community) & is.null(vertex)) {
    p <- p + geom_node_point(aes_string(color = "comm_id"), size = vertex_size) +
      scale_color_manual(limits = names(pal_community), values = pal_community)
  } else if(is.null(community) & !is.null(vertex)) {
    p <- p + geom_node_point(aes_string(color = vertex), size = vertex_size) +
      scale_color_manual(limits = names(vertex_pal), values = vertex_pal)
  } else if(is.null(community) & is.null(vertex)) {
    p <- p + geom_node_point(color = "gray65", size = vertex_size) 
  }
  if(vertex_label) {
    p <- p +
      geom_node_text(aes_string(label = "label"), repel = T)
  }
  if(!is.null(file_out)) {
    jpeg(file_out, res = res, height = height, width = width, units = units)
    print(p)
    dev.off()
    return( ct.g)
  } else {
    return(list(plot = p, graph = ct.g))
  }
  
}
