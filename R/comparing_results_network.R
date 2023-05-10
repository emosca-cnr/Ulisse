#' Visual representation of multiple cross-talk analyses.
#' @description The function is used to merge the results obtained by multiple cross-talk analyses by creating a union network and highlighting the
#' shared edges and vertices
#' @details The function takes as an input a list of data.frames resulting from multiple `gs_cross_talk()` analyses
#' and produce a union network that will be plotted with the edges colored by the number or which of results that have that cross-talk 
#' @param res_list list of resulting data.frame from `gs_cross_talk()`
#' @param vertex_number logical, if the vertices should be colored by the number of results in which are present
#' @param vertex_number_pal vector of two colors to be used for gradient coloring the vertices when `vertex_number = TRUE`. 
#'  If `NULL` the gradient is build from "red" to "blue".
#' @param vertex_number_adj value used to set transparency to vertex colors in `vertex_number_pal`
#' @param vertex_size size of the vertices
#' @param vertex can be either `name` value to color the vertices according to their name, or a list of two vector, where the first is 
#'  the name of the attribute and the second is a named vector with the attribute per vertices 
#'  to be used for color (only discrete value), or `NULL` if the vertices should be not colored (thus "grey65" is passed as default color). 
#'  If `vertex_number = TRUE`, `vertex` is represented as cells of voronoi tessellation (see `ggraph::geom_node_voronoi()`)
#' @param vertex_pal named vector with the colors to be used for each unique `vertex` element. If `NULL` pals::alphabet2() palette is used
#' @param voronoi_radius,voronois_alpha: parameters passed to `geom_node_voronoi()`. If both `vertex` and `vertex_number` arguments are enabled, then the 
#'  first is represented as cells of voronoi tessellation (see `ggraph::geom_node_voronoi()`). In this case, these argument are used to control 
#'  the radius of the cells and the transparency, respectively.
#' @param edge_color_by=c("number","which") if the edges should be colored by number of results that share that edge ("number") or which results contain it ("which")
#' @param edge_width logical, if the edge width should be proportional to the number of results that share that edge. 
#' @param edge_pal palette used to color the edges according to `edge_color_by`. If `edge_color_by="number` than the palette should be the two color used for
#'  gradient coloring them; if `NULL` "red" to "blue" is used. If `edge_color_by="which` than the palette should be a (named) vector with enought colors for each 
#'  results name combination
#' @param file_out name used to save the plot in jpeg format plot. If `NULL` the functions returns also the plot object
#' @param width,height,res,units graphical value of `jpeg()` function
#' @param ... further graphical parameters to be passed to `ggraph()` function
#' @return If `file_out` is null the function returns the plot and the igraph object used for plotting. 
#'  Otherwise, only the igraph object is returned and the plot is saved through to `file_out` 
#' @importFrom grDevices rainbow dev.off jpeg adjustcolor
#' @import igraph
#' @import ggraph
#' @export




comparing_results_network <- function(res_list, file_out, vertex_number = F, vertex_number_pal = NULL, vertex_number_adj = 0.7, vertex_size = 5,
                                      vertex = NULL, vertex_pal =NULL, voronoi_radius = 0.8, voronoi_alpha = 0.3, vertex_label = TRUE,
                                      edge_color_by = "number", edge_width = TRUE,
                                      edge_pal = NULL, edge_adj_col = 0.7, 
                                      width = 200, height = 200, res = 300, units = "mm",  ...) {
  
  res_list <- lapply(res_list, function(x) {
    x$sample <- 1
    x <- graph_from_data_frame(x, directed = F)
    V(x)$sample <- 1
    return(x)
  }
  )
  n <- names(res_list)
  res_list <- lapply(1:length(res_list), function(x) {
    E(res_list[[x]])$s_name <- n[x]
    return(res_list[[x]])
  })
  names(res_list) <- n
  
  union.g <- igraph::union(res_list[[1]], res_list[2:length(res_list)])
  attr.union.g <- data.frame(edge.attributes(union.g), stringsAsFactors = F)
  
  idx <- grep("sample", colnames(attr.union.g))
  attr.union.g$nsample <- rowSums(attr.union.g[, idx], na.rm = T)
  
  E(union.g)$edge_number <- attr.union.g$nsample
  
  if(edge_color_by == "number") {
    edge_color_by = "edge_number"
  } else if(edge_color_by == "which") {
   idx <- grep("s_name", colnames(attr.union.g))
   attr.union.g$res <- apply( attr.union.g[ , idx ] , 1 , function(x) paste(x[!is.na(x)], collapse = ";"))
   E(union.g)$res <- attr.union.g$res
  }
  
  if(edge_width) {
    edge_width = "edge_number"
  }
  #building defaults and value to be plotted---------------
  #filename
  # if(is.null(file_out)) {
  #   file_out <- "union_ct_network.jpeg"
  # }
  #vertex color-------------------------------------
  if(vertex_number) {
    Vattr.union.g <- data.frame(vertex.attributes(union.g), stringsAsFactors = F)
    
    idx <- grep("sample", colnames(Vattr.union.g))
    Vattr.union.g$nsample <- rowSums(Vattr.union.g[, idx], na.rm = T)
    
    V(union.g)$vertex_number <- Vattr.union.g$nsample
    
    if(is.null(vertex_number_pal)) {
      vertex_number_pal <- c("blue", "red")
    }
  }
  
  #other possible annotations-- if both number and annotations --> voronoi
  if(!is.null(vertex)) {
    if(length(vertex) == 2) {
      vertex_attr(union.g, name = vertex[[1]]) <- vertex[[2]][V(union.g)$name]
      vertex <- vertex[[1]]
    } 
  }
  if(!is.null(vertex) & is.null(vertex_pal)) {
    v <- unique(vertex_attr(union.g, name = vertex))
    vertex_pal <- setNames(pals::alphabet2(length(v)), v)
  }
  
  #labels
  if(is.logical(vertex_label)) {
    if(vertex_label) {
      V(union.g)$label <- V(union.g)$name
    }
  } else {
    V(union.g)$label <- vertex_label[V(union.g)$name]
    vertex_label <- TRUE
  }
  #edges color-----------------------------------
  
  if(is.null(edge_pal)) {
    if(edge_color_by == "edge_number") {
      edge_pal <- c("blue", "red")
    } else if (edge_color_by == "which") {
      ed <- unique(E(union.g)$res)
      edge_pal <- setNames(pals::alphabet2(length(ed)), ed)
    }
    
  } else {
    if(edge_color_by == "which") {
      ed <- unique(E(union.g)$res)
      n_ed <- names(edge_pal)
      idx <- sum(n_ed %in% ed)
      if(length(edge_pal) >= length(ed) & idx == length(ed)) {
        edge_pal <- setNames(edge_pal, ed)
      } else {
        print.warnings("Not enought color provided as edge_pal palette. Using pals::alphabet2() instead")
        edge_pal <- setNames(pals::alphabet2(length(ed)), ed)
      }
      
    }
  }
  
  
  #plotting----------------------------
  
  p <- ggraph(union.g, ...) +
    theme_graph() 
  if(vertex_number & !is.null(vertex)) {
    p <- p + geom_node_voronoi(aes_string(fill = vertex), max.radius = voronoi_radius, colour = 'white', alpha = voronoi_alpha) +
      scale_fill_manual(limits = names(vertex_pal), values = vertex_pal) 
  }
  
  if(!is.null(edge_width) & edge_color_by == "edge_number") {
    p <- p + geom_edge_link(aes_string(edge_colour = edge_color_by, edge_width = edge_width), alpha = edge_adj_col) +
      scale_edge_color_gradient(low = edge_pal[1], high = edge_pal[2])
    
  } else if(!is.null(edge_width) & edge_color_by == "which") {
    p <- p + geom_edge_link(aes_string(edge_colour = "res", edge_width = edge_width), alpha = edge_adj_col) +
      scale_edge_color_manual(limits = names(edge_pal), values = edge_pal)
    
  } else if(is.null(edge_width) & edge_color_by == "edge_number") {
    p <- p + geom_edge_link(aes_string(edge_colour = edge_color_by), alpha = edge_adj_col) +
      scale_edge_color_gradient(low = edge_pal[1], high = edge_pal[2])
    
  } else if (is.null(edge_width) & edge_color_by == "which") {
    p <- p + geom_edge_link(aes_string(edge_colour = res), alpha = edge_adj_col) +
      scale_edge_color_manual(limits = names(edge_pal), values = edge_pal)
  }
  
  if(vertex_number & !is.null(vertex)) {
    p <- p + geom_node_point(aes(color = vertex_number), size = vertex_size) +
      scale_color_gradient(low = vertex_number_pal[1], high = vertex_number_pal[2])
    
  } else if(vertex_number  & is.null(vertex)) {
    p <- p + geom_node_point(aes_string(color = vertex_number), size = vertex_size) +
      scale_color_gradient(low = vertex_number_pal[1], high = vertex_number_pal[2])
  } else if(!vertex_number & !is.null(vertex)) {
    p <- p + geom_node_point(aes_string(color = vertex), size = vertex_size) +
      scale_color_manual(limits = names(vertex_pal), values = vertex_pal)
  } else if(is.null(community) & is.null(vertex)) {
    p <- p + geom_node_point(color = "gray65", size = vertex_size) 
  }
  if(vertex_label) {
    p <- p +
      geom_node_text(aes(label = label), repel = T)
  }
  
  if(!is.null(file_out)) {
    jpeg(file_out, res = res, height = height, width = width, units = units)
    print(p)
    dev.off()
    return( union.g)
  } else {
    return(list(plot = p, graph = union.g))
  }
  
  
}


