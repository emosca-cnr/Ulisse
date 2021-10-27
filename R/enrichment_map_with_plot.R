#' Enrichment map
#' @param x named vector of pathway scores
#' @param gs_list gene set list
#' @param method 'overlap' or 'jaccard'
#' @param coeff threshold for the similarity score between two gene sets
#' @param all_gs TRUE/FALSE
#' @param weight.within value to weight the attraction between two vertices of the same community
#' @param weight.between value to weight the attraction between two vertices of two disntict communities
#' @param file output file
#' @param comm_method one among 'auto', 'fastgreedy' and 'multilev'
#' @param vertex.label.cex vertex label size
#' @param vertex.label.dist vertex label distance from veritex
#' @param vertex.label.degree vertex label orientation
#' @param vertex.label.font vertex label font
#' @param gs_list_size gene set list size
#' @param set_sim_df optional data frame with three columns 'set1', 'set2' and 'sim'
#' @param plot_flag TRUE/FALSE
#' @param vertex.color vertex color
#' @param img.width image width
#' @param img.height image height
#' @param img.res image resolution
#' @param min_comm_size min community size
#' @param vertex.size.min min vertex size
#' @param vertex.size.max max vertex size
#' @param edge.wd.min min edge widht
#' @param edge.wd.max max edge width
#' @param mark.groups TRUE/FALSE
#' @param vertex.shape vertex shape
#' @param vertex.pie vertex pie
#' @param ... additional parameters of plot.igraph
#' @return list of two data frames containing pathway network and vertex similarity
#' @import igraph
#' @importFrom  RColorBrewer brewer.pal
#' @export
#'
enrichment_map_with_plot <- function(x, gs_list, method=c('overlap', 'jaccard'), coeff=c(0.5, 0.1), all_gs=TRUE, weight.within = 4, weight.between = 1, file='en_map.jpg', comm_method=c('auto', 'fastgreedy', 'multilev'), vertex.label.cex=1, vertex.label.dist=0.5, vertex.label.degree=pi/2, vertex.label.font=2, gs_list_size=NULL, set_sim_df=NULL, plot_flag=TRUE, vertex.color=NULL, img.width=200, img.height=200, img.res=300, min_comm_size=2, vertex.size.min=2, vertex.size.max=5, edge.wd.min=1, edge.wd.max=4, mark.groups=T, vertex.shape='circle', vertex.pie=NULL, ...){


  method <- match.arg(method)
  comm_method <-  match.arg(comm_method)
  if(method == "overlap"){
    coeff <- coeff[1]
  }
  if(method == 'jaccard'){
    corff <- coeff[2]
  }

  if(is.null(gs_list_size)){
    gs_list_size <- unlist(lapply(gs_list, length))
    gs_list_subset <- gs_list[names(gs_list) %in% names(x)]
  }

  if(is.null(set_sim_df)){
    set_sim_df <- data.frame(t(combn(names(x), 2)), stringsAsFactors = FALSE)
    colnames(set_sim_df)[1:2] <- c('set1', 'set2')
    for(j in 1:nrow(set_sim_df)){
      if(j %% 1000 ==0 ){
        cat(j, '/', nrow(set_sim_df), '\n')
      }
      if(method == 'overlap'){
        set_sim_df$sim[j] <- calc_set_similarity(gs_list_subset[names(gs_list_subset) == set_sim_df$set1[j]][[1]], gs_list_subset[names(gs_list_subset) == set_sim_df$set2[j]][[1]], method = 'overlap')
      }
      if(method == 'jaccard'){
        set_sim_df$sim[j] <- calc_set_similarity(gs_list_subset[names(gs_list_subset) == set_sim_df$set1[j]][[1]], gs_list_subset[names(gs_list_subset) == set_sim_df$set2[j]][[1]], method='jaccard')
      }
    }
  }else{
    gs_list_subset <- unique(c(set_sim_df[, 1], set_sim_df[, 2]))
    gs_list_subset <- gs_list_subset[ gs_list_subset %in% names(x)]
  }


  #enrichment map: overlap
  path_mod <- igraph::simplify(igraph::graph.data.frame(set_sim_df[set_sim_df$sim >= coeff, ], directed=FALSE), edge.attr.comb = 'mean')
  if(length(igraph::V(path_mod)) < length(x) & all_gs){
    path_mod <- igraph::add_vertices(path_mod, length(which(!(names(x) %in% igraph::V(path_mod)$name))), attr = list(name=names(x)[which(!(names(x) %in% igraph::V(path_mod)$name))])) #reintroduce vertexes excluded
  }

  #community detection
  path_mod_comm <- find_communities(path_mod, e.weights = igraph::E(path_mod)$sim, methods = c('fastgreedy', 'multilev'))
  print(path_mod_comm$info)
  if(comm_method == 'auto'){
    path_mod_comm <- path_mod_comm$comm[names(path_mod_comm$comm) == path_mod_comm$info$algorithm[which.max(path_mod_comm$info$modularity)]][[1]]
  }
  if(comm_method == 'fastgreedy'){
    path_mod_comm <- path_mod_comm$comm$fastgreedy
  }
  if(comm_method == 'multilev'){
    path_mod_comm <- path_mod_comm$comm$multilev
  }


  out <- data.frame(id=names(x), score=x, module=0, size=1, stringsAsFactors = FALSE)
  out$module[out$id %in% path_mod_comm$names]<- path_mod_comm$membership[match(out$id[out$id %in% path_mod_comm$names], path_mod_comm$names)]
  out$size[out$id %in% path_mod_comm$names] <- table(path_mod_comm$membership)[out$module[out$id %in% path_mod_comm$names]]

  if(min_comm_size > 1){
    out <- out[out$size >=min_comm_size, ]
    path_mod <- igraph::induced_subgraph(path_mod, igraph::V(path_mod)$name[igraph::V(path_mod)$name %in% out$id[out$size >=min_comm_size]])
    out <- out[match(igraph::V(path_mod)$name, out$id), ]
  }

  #print pathway network - option 1
  if(plot_flag){

    #layout
    L <- igraph::layout_with_fr(path_mod, weights = edge_weights(list(membership=out$module), path_mod, weight.within = weight.within, weight.between = weight.between))

    #vertex color
    if(vertex.shape != 'pie'){
      if(is.null(vertex.color)){
        igraph::V(path_mod)$color <- RColorBrewer::brewer.pal(9, "Greens")[round(linear_map(x[match(igraph::V(path_mod)$name, names(x))], 1, 9))]
      }else{
        igraph::V(path_mod)$color <- vertex.color[match(igraph::V(path_mod)$name, names(vertex.color))]
      }
    }
    #vertex shape
    if(length(vertex.shape) == length(V(path_mod))){
      igraph::V(path_mod)$shape <- vertex.shape[match(igraph::V(path_mod)$name, names(vertex.shape))]
    }
    if(vertex.shape == 'pie'){
      igraph::V(path_mod)$shape <- 'pie'
      igraph::V(path_mod)$pie <- vertex.pie[match(igraph::V(path_mod)$name, names(vertex.pie))]
    }
    temp <- cbind(L, comm=out$module, id=1:nrow(L))
    temp <- temp[order(temp[, 3], temp[,2]), ]
    temp <- temp[!duplicated(temp[, 3]), ]
    igraph::V(path_mod)$label <- ''
    igraph::V(path_mod)$size <- round(linear_map(gs_list_size[match(igraph::V(path_mod)$name, names(gs_list_size))], vertex.size.min, vertex.size.max))
    igraph::E(path_mod)$width <- round(linear_map(igraph::E(path_mod)$sim, edge.wd.min, edge.wd.max))

    if(mark.groups){
      mark.groups <- split(igraph::V(path_mod)$name, out$module)
    }else{
      mark.groups <- NULL
    }

    #without labels
    jpeg(file, units='mm', width = img.width, height = img.height, res=img.res)
    par(mar=c(2, 0, 0, 0), oma=c(0, 0, 0, 0))
    igraph::plot.igraph(path_mod, mark.groups = mark.groups, mark.expand = 2, layout=L, ...)

    if(is.null(vertex.color) & vertex.shape != 'pie'){
      legend_text <- c(round(sort(x)[1], 2), rep(NA, 3), round(median(x), 2), rep(NA, 3), round(sort(x)[length(x)], 2))
      legend("bottomright", legend = legend_text, pch = 22, pt.bg = RColorBrewer::brewer.pal(9, "Greens"), col = "black", xpd = TRUE, bty = "n", cex = 0.8)
    }

    dev.off()

    #with LABELS
    igraph::V(path_mod)$label[temp[, 4]] <- temp[, 3]

    jpeg(paste(gsub(".jpg", "_labels", file), ".jpg", sep=''), units='mm', width = img.width, height = img.height, res=img.res)
    par(mar=c(2, 0, 0, 0), oma=c(0, 0, 0, 0))
    igraph::plot.igraph(path_mod, mark.groups = mark.groups, layout=L, vertex.label.cex=vertex.label.cex, vertex.label.dist=vertex.label.dist, vertex.label.degree=vertex.label.degree, vertex.label.font=vertex.label.font, vertex.pie=V(path_mod)$pie, ...)

    if(is.null(vertex.color) & vertex.shape != 'pie'){
      legend_text <- c(round(sort(x)[1], 2), rep(NA, 3), round(median(x), 2), rep(NA, 3), round(sort(x)[length(x)], 2))
      legend("bottomright", legend = legend_text, pch = 22, pt.bg = RColorBrewer::brewer.pal(9, "Greens"), col = "black", xpd = TRUE, bty = "n", cex = 0.8)
    }

    dev.off()
  }
  return(list(network=out, sim_coeff=set_sim_df))

}
