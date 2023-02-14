#' Function for plotting results from PCT or CCC analysis as a network
#' @description The function elaborates the results of `pathway_cross_talk()` or `cluster_communication()` functions to obtain a graphical 
#' representation of the cross-talk network
#' @details The functions uses `pathway_cross_talk()` or `cluster_communication()` output to build a CT network. Then, the network is plotted with
#' the edges colored by a value of interest. The user may decide if to plot all the edges, color the vertices by 
#' their community, scale the dimension of the edges over a variable of interest.
#' @param ct output of `pathway_cross_talk()` or `cluster_communication()` function
#' @param all logical, if the function have to plot all the edges or only the significant ones. If `TRUE`, the function uses `p_val, FDR` and `ct_val` to identify the 
#' significant ones, that are plotted with a solid line, and the not significant with a dashed line. If `FALSE`, the functions uses the same filtering parameters to 
#' filter the data and plot only the CT passing the filtering threshold as a solid line. If `NULL` the function does not 
#' filter the data and plot all the edges as a solid line
#' @param p_val,FDR,ct_val filtering values. If one of these is set to `NULL` the function ignores it
#' @param community logical or an object resulting from igraph community calculation. If `TRUE`, the function calculates
#' the communities by using fastgreedy algorithm, if `FALSE` the function does not consider communities. Otherwise, if a
#' community object is provided, the function uses it and color the vertices accordingly
#' @param pal_community color palette to be used to color the vertices according to the communities. If not provided,
#' the function uses `rainbow()` palette
#' @param edge_col_by name of the column to be used to color the edges. If `NULL` the edges will be not colored
#' @param edge_pal color to be used to create the gradient to color the edges. If not provided, the function uses "blue" 
#' to "red"
#' @param edge_value vector of the two values to be used to create the gradient to color the edges. If not provided, 
#' the function uses minimum and maximum
#' @param edge_width logical, if the edge widths should be proportional to `edge_col_by` values. If `TRUE` the `edge_col_by` values are scaled between 1 to 9 to set the
#' edge width accordingly
#' @param edge_adj_col value used to adjust color transparency of the edges
#' @param file_out name of the jpeg file produced. If not provided, the plot will be automatically named `ct_network.jpeg`
#' @param width,height,res,units graphical value of `jpeg()` function
#' @param layout layout to be used in igraph plotting. If not provided, the function calculates it
#' @param legend_edge_col logical, if to plot the legend of the color of the edges
#' @param legend_filt logical, if to plot legend of the different lines used when `all = T`
#' @param legend_community logical, if to plot the legend of the community color used
#' @param x_filt,y_filt,x_comm,y_comm, x and y coordinates for `legend_filt`, `legend_community`, 
#' positioning in the plot. See `legend()` function for more details
#' @param xl,yb,xr,yt coordinates for the positioning of the gradient legend for `legend_edge_col`. See `color.legend()` 
#' and `rect()` functions for more details
#' @param ... further graphical parameters to be passed to `igraph.plot()` function
#' @return the function produce the plot saved with the name passed to `file_out` and also returns the layout,
#' the igraph object used for plotting and the communities (if calculated)
#' @importFrom graphics legend
#' @importFrom grDevices rainbow dev.off jpeg adjustcolor
#' @importFrom plotrix color.legend rescale
#' @importFrom circlize colorRamp2
#' @import igraph
#' @export


plot_network_CT <- function(ct, all = FALSE, p_val, FDR, ct_val, community, pal_community =NULL, 
                             edge_col_by = "pct", edge_pal=NULL, edge_value =NULL, edge_width = T, edge_adj_col = 0.7,
                             file_out =NULL,  width = 200, height = 200, res = 300, units = "mm", 
                             layout = NULL, legend_edge_col=T, legend_filt = T, legend_community = T, 
                             x_filt = -1.5, y_filt = 0, 
                             x_comm = -1.5, y_comm = -0.5, 
                             xl = -1.45, yb = 0.2, xr = -1.3, yt = 0.5, ...) {
  
  #building defaults and values to be plotted--------------
  #filename
  if(is.null(file_out)) {
    file_out <- "ct_network.jpeg"
  }
  ct_column <- grep("pct|ccc_score", colnames(ct))
  ct_column <- colnames(ct)[ct_column]
  #filtering
  if(is.null(all)) {
    legend_filt <- F
    ct$filt <- 1
  } else if(!all) {
    if(is.null(p_val)) {
      p_val <- max(ct$p_value_link, na.rm = T)
    }
    if(is.null(FDR)) {
      FDR <- max(ct$FDR_link, na.rm = T)
    }
    if(is.null(ct_val)) {
      ct_val <- min(ct[,ct_column], na.rm = T)
    }
    ct <- ct[ct$p_value_link <= p_val & ct$FDR_link <= FDR & ct[,ct_column] >= ct_val,]
    ct$filt <- 1
    legend_filt <- F
  } else if(all) {
    if(is.null(p_val)) {
      p_val <- max(ct$p_value_link, na.rm = T)
    }
    if(is.null(FDR)) {
      FDR <- max(ct$FDR_link, na.rm = T)
    }
    if(is.null(ct_val)) {
      ct_val <- min(ct[,ct_column], na.rm = T)
    }
    ct$filt <- 3
    ct$filt[ct$p_value_link <= p_val & ct$FDR_link <= FDR & ct[,ct_column] >= ct_val] <- 1
  }
  #network
  ct.g <- graph_from_data_frame(ct, directed = F)
  #community
  if(is.logical(community)) {
    if(!community) {
      community <- NULL
    } else {
      community <- fastgreedy.community(ct.g, weights = edge_attr(ct.g)[[ct_column]])
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
  #edges
  if(is.null(edge_pal)) {
    edge_pal <- c("blue", "red")
  }
  
  if(!is.null(edge_col_by)) {
    edge_att <- edge_attr(ct.g, name = edge_col_by)
    if(is.null(edge_value)) {
      
      edge_value <- c(min(edge_att), max(edge_att))
    }
  }
  if(edge_width) {
    if(!is.null(edge_col_by)) {
      edge_width <- rescale(edge_att, c(1,9))
    } else {
      edge_att <- edge_attr(ct.g, name = ct_column)
      edge_width <- rescale(edge_att, c(1,9))
    }
  } else {
    edge_width <- NULL
    
  }
  
  
  if(is.null(layout)) {
    lo_cl <- layout_with_fr(ct.g, weights = edge_att)
  } else {
    lo_cl <- layout
  }
  
  #plotting-----------------------------
  edge_cl <- colorRamp2(edge_value, edge_pal)
  E(ct.g)$color <- adjustcolor(edge_cl(edge_att), edge_adj_col)
  
  
  if(is.null(community)) {
    if(is.null(edge_width)) {
      if(is.null(edge_col_by)) {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, edge.lty = E(ct.g)$filt, layout=lo_cl, ...)
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not Significant"), box.col = NA)
        }
        
        
        dev.off() 
      } else {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, edge.col = E(ct.g)$color, edge.lty = E(ct.g)$filt, layout=lo_cl, ...)
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not significant"), box.col = NA)
        }
        if(legend_edge_col) {
          color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                       legend = c(round(edge_value[1], y = 3), round((edge_value[1] + edge_value[2])/2, y = 3), round(edge_value[2], y = 3)), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        
        dev.off()
      }
    } else {
      if(is.null(edge_col_by)) {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, edge.lty = E(ct.g)$filt,
             edge.width = edge_width, layout=lo_cl, ...)
        
        
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not Significant"), box.col = NA)
        }
        
        dev.off()
      } else {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, edge.col = E(ct.g)$color, edge.lty = E(ct.g)$filt,
             edge.width = edge_width, layout=lo_cl, ...)
        
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not Significant"), box.col = NA)
        }
        
        if(legend_edge_col) {
          color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                       legend = c(round(edge_value[1], y = 3), round((edge_value[1] + edge_value[2])/2, y = 3), round(edge_value[2], y = 3)), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        dev.off()
      }
    }
  } else {
    if(is.null(edge_breaks)) {
      if(is.null(edge_col_by)) {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, vertex.color=pal_community[as.character(comm_ct_mem)],
             edge.lty = E(ct.g)$filt, layout=lo_cl, ...)
        if(legend_community) {
          legend(x = x_comm, y = y_comm, fill = pal_community, legend = names(pal_community))
        }
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not Significant"), box.col = NA)
        }
        
        dev.off()
      } else {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, edge.col = E(ct.g)$color, vertex.color=pal_community[as.character(comm_ct_mem)],
             edge.lty = E(ct.g)$filt, layout=lo_cl, ...)
        if(legend_community) {
          legend(x = x_comm, y = y_comm, fill = pal_community, legend = names(pal_community))
        }
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not Significant"), box.col = NA)
        }
        if(legend_edge_col) {
          color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                       legend = c(round(edge_value[1], y = 3), round((edge_value[1] + edge_value[2])/2, y = 3), round(edge_value[2], y = 3)), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        
        dev.off()
      }
    } else {
      if(is.null(edge_col_by)) {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, edge.lty = E(ct.g)$filt, vertex.color=pal_community[as.character(comm_ct_mem)],
             edge.width = edge_width, layout=lo_cl, ...)
        if(legend_community) {
          legend(x = x_comm, y = y_comm, fill = pal_community, legend = names(pal_community))
        }
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not Significant"), box.col = NA)
        }
        
        dev.off()
      } else {
        jpeg(file_out, width = width, height = height, res=res, 
             units=units)
        plot(ct.g, edge.col = E(ct.g)$color, vertex.color=pal_community[as.character(comm_ct_mem)],
             edge.width = edge_width, edge.lty = E(ct.g)$filt, layout=lo_cl, ...)
        if(legend_community) {
          legend(x = x_comm, y = y_comm, fill = pal_community, legend = names(pal_community))
        }
        
        if(legend_filt) {
          legend(x= x_filt, y = y_filt, col = "black", lty = c(1,3), lwd = 3,
                 legend = c("Significant", "Not Significant"), box.col = NA)
        }
        if(legend_edge_col) {
          color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                       legend = c(round(edge_value[1], y = 3), round((edge_value[1] + edge_value[2])/2, y = 3), round(edge_value[2], y = 3)), 
                       rect.col = edge_cl(seq(edge_value[1], edge_value[2], length.out = 10 )), 
                       gradient="y", align = "rb")
        }
        
        dev.off()
      }
    }
  }
  
  #returning objects
  if(!is.null(community)) {
    return(list(layout = lo_cl, graph = ct.g, community = community))
  } else {
    return(list(layout = lo_cl, graph = ct.g))
  }
  
  
  
  
}
