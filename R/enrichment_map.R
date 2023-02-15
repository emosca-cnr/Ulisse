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
#' @import igraph pals
#' @importFrom  RColorBrewer brewer.pal
#' @export
#'
enrichment_map <- function(x, gs_list, method=c('overlap', 'jaccard'), coeff=NULL, all_gs=TRUE, weight.within = 4, weight.between = 1, file.prefix='en_map', comm_method='fastgreedy', vertex.label.cex=1, vertex.label.dist=0.5, vertex.label.degree=pi/2, vertex.label.font=2, gs_list_size=NULL, set_sim_df=NULL, plot_flag=TRUE, vertex.color=NULL, img.width=200, img.height=200, img.res=300, min_comm_size=2, vertex.size.min=2, vertex.size.max=5, edge.wd.min=1, edge.wd.max=4, mark.groups=T, vertex.shape='circle', vertex.pie=NULL, vertex.color.pal=NULL, n_width=3, n_size=3, L=NULL, score.decreasing=TRUE, n_name=1, legend.cex=0.7, ...){
  
  
  method <- match.arg(method)
  cat("method:", method, "\n")
  if(is.null(coeff)){
    if(method == "overlap"){
      coeff <- 0.5
    }
    if(method == 'jaccard'){
      coeff <- 0.1
    }
  }
  cat("coeff:", coeff, "\n")
  cat("comm_method:", comm_method, "\n")
  
  if(is.null(gs_list_size)){
    gs_list <- gs_list[names(gs_list) %in% names(x)]
    gs_list_size <-lengths(gs_list)
  }
  
  if(is.null(set_sim_df)){
    cat("Calculating gene set similarities...")
    set_sim_df <- data.frame(t(utils::combn(names(x), 2)), stringsAsFactors = FALSE)
    colnames(set_sim_df)[1:2] <- c('set1', 'set2')
    for(j in 1:nrow(set_sim_df)){
      if(j %% 1000 ==0 ){
        cat(j, '/', nrow(set_sim_df), '\n')
      }
      set_sim_df$sim[j] <- calc_set_similarity(gs_list[names(gs_list) == set_sim_df$set1[j]][[1]], gs_list[names(gs_list) == set_sim_df$set2[j]][[1]], method=method)
    }
    cat("done.\n")
  }else{
    gs_list <- unique(c(set_sim_df[, 1], set_sim_df[, 2]))
    gs_list <- gs_list[ gs_list %in% names(x)]
  }
  
  
  cat("Summary of gene set similarities (", nrow(set_sim_df),"):\n")
  print(summary(set_sim_df$sim))
  
  #enrichment map: overlap
  path_mod <- igraph::simplify(igraph::graph.data.frame(set_sim_df[set_sim_df$sim >= coeff, ], directed=FALSE), edge.attr.comb = 'mean')
  if(length(igraph::V(path_mod)) < length(x) & all_gs){
    path_mod <- igraph::add_vertices(path_mod, length(which(!(names(x) %in% igraph::V(path_mod)$name))), attr = list(name=names(x)[which(!(names(x) %in% igraph::V(path_mod)$name))])) #reintroduce vertexes excluded
  }
  cat("Enrichment Map: V=", length(V(path_mod)), ", E=", length(E(path_mod)), "\n")
  
  cat("Community detection...\n")
  path_mod_comm <- find_communities(path_mod, e.weights = igraph::E(path_mod)$sim, methods = comm_method)
  print(path_mod_comm$info)
  
  path_mod_comm <- path_mod_comm$comm[[comm_method]]
  
  V(path_mod)$comm_id <- path_mod_comm$membership
  # out <- data.frame(id=names(x), score=x, module=0, size=1, stringsAsFactors = FALSE)
  # out$module[out$id %in% path_mod_comm$names]<- path_mod_comm$membership[match(out$id[out$id %in% path_mod_comm$names], path_mod_comm$names)]
  # out$size[out$id %in% path_mod_comm$names] <- table(path_mod_comm$membership)[out$module[out$id %in% path_mod_comm$names]]
  
  if(min_comm_size > 1){
    cat("Excluding communities smaller than ", min_comm_size, "\n")
    keep_comm <- sizes(path_mod_comm)
    keep_comm <- names(keep_comm)[keep_comm >= min_comm_size]
    cat("Kept:", length(keep_comm), "\n")
    path_mod <- igraph::induced_subgraph(path_mod, igraph::V(path_mod)$name[igraph::V(path_mod)$comm_id %in% keep_comm])
    cat("Enrichment Map: V=", length(V(path_mod)), ", E=", length(E(path_mod)), "\n")
  }
  
  V(path_mod)$score <- x[match(V(path_mod)$name, names(x))]
  
  ### genes of each pathway community
  path_comm_genes <- split(V(path_mod)$name, V(path_mod)$comm_id)
  path_comm_genes <- lapply(path_comm_genes, function(pathway_vector) sort(unique(unlist(gs_list[names(gs_list) %in% pathway_vector]))))
  
  #print pathway network 
  if(plot_flag){
    
    if(is.null(vertex.color.pal)){
      vertex.color.pal <- pals::brewer.greens(n = 5)
      vertex.color.factor <- ggplot2::cut_interval(V(path_mod)$score, n=5)
    }
    
    #color
    mark.col <- pals::polychrome(length(table(V(path_mod)$comm_id)))
    
    #layout
    if(is.null(L)){
      L <- igraph::layout_with_fr(path_mod, weights = edge_weights(list(membership=V(path_mod)$comm_id), path_mod, weight.within = weight.within, weight.between = weight.between))
    }
    
    #vertex color
    if(vertex.shape != 'pie'){
      if(is.null(vertex.color)){
        igraph::V(path_mod)$color <- vertex.color.pal[as.numeric(vertex.color.factor)]
      }else{
        igraph::V(path_mod)$color <- vertex.color
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
    
    ### community labels
    network_df <- data.frame(L, name=V(path_mod)$name, comm_id=V(path_mod)$comm_id, score=V(path_mod)$score, show.name=FALSE, stringsAsFactors = F) ##all
    if(score.decreasing){
      temp <- tapply(setNames(-network_df$score, network_df$name), network_df$comm_id, rank, ties.method="min")
    }else{
      temp <- tapply(setNames(network_df$score, network_df$name), network_df$comm_id, rank, ties.method="min")
    }
    temp <- data.frame(rank=unlist(temp), name=unlist(lapply(temp, names)), stringsAsFactors = F)
    network_df <- merge(network_df, temp, by="name")
    rm(temp)
    network_df$show.name[network_df$rank <= n_name] <- TRUE
    
    #community label"
    comm_top <- data.frame(show.comm=unlist(tapply(network_df$X2, network_df$comm_id, max)))
    network_df <- merge(network_df, comm_top, by.x="comm_id", by.y=0, all=T)
    network_df$show.comm <- network_df$X2 == network_df$show.comm
    V(path_mod)$label <- network_df$show.comm[match(V(path_mod)$name, network_df$name)]
    V(path_mod)$label[V(path_mod)$label] <- V(path_mod)$comm_id[V(path_mod)$label]
    V(path_mod)$label[V(path_mod)$label == FALSE] <- ""
    
    network_df$color <- mark.col[as.numeric(as.factor(network_df$comm_id))]
    
    size.pal <- seq(from=vertex.size.min, to=vertex.size.max, length.out=n_size)
    size.factor <- ggplot2::cut_interval(gs_list_size[match(igraph::V(path_mod)$name, names(gs_list_size))], n = n_size)
    igraph::V(path_mod)$size <- size.pal[as.numeric(size.factor)]
    
    width.pal <- seq(from=edge.wd.min, to=edge.wd.max, length.out=n_width)
    width.factor <- ggplot2::cut_interval(igraph::E(path_mod)$sim, n = n_width)
    igraph::E(path_mod)$width <- width.pal[as.numeric(width.factor)]
    
    ### mark groups
    if(mark.groups){
      mark.groups <- split(igraph::V(path_mod)$name, V(path_mod)$comm_id)
    }else{
      mark.groups <- NULL
    }
    
    #without labels
    grDevices::jpeg(paste0(file.prefix, ".jpg"), units='mm', width = img.width, height = img.height, res=img.res)
    
    layout(matrix(c(2, 2, 1, 3), nrow = 2, byrow = T), heights = c(0.10, 0.90), widths = c(.8, .2))
    
    par(mar=c(.1, .1, .1, .1), oma=c(0, 0, 0, 0))
    igraph::plot.igraph(path_mod, mark.groups = mark.groups, layout=L, mark.col=adjustcolor(mark.col, 0.2), mark.border=NA, vertex.label.degree=-pi/2, vertex.label.font=2, vertex.label.cex=2, vertex.label.dist=1, vertex.label.color="black", ...)
    
    ### module annotations
    par(mar=c(0.1, 0.1, 0.1, 0.1))
    plot.new()
    #legend("topleft", legend = paste(network_df$comm_id[network_df$show.name], network_df$name[network_df$show.name]), cex=0.7, text.col = network_df$color[network_df$show.name], bty = "n", xpd=T, text.font = 2, ncol=2)
    legend("topleft", legend = paste(network_df$comm_id[network_df$show.name], network_df$name[network_df$show.name]), cex=legend.cex, bty = "n", xpd=T, text.font = 2, ncol=2)
    
    if(is.null(vertex.color) & vertex.shape != 'pie'){
      plot.new()
      graphics::legend("topright", legend = levels(size.factor), pch = 1, pt.cex =size.pal, xpd = TRUE, bty = "n", cex = 0.8)
      
      graphics::legend("right", legend = levels(width.factor), lwd=width.pal, xpd = TRUE, bty = "n", cex = 0.8)
      
      graphics::legend("bottomright", legend = levels(vertex.color.factor), pch = 21, pt.bg = vertex.color.pal, col = "black", xpd = TRUE, bty = "n", cex = 0.8)
    }
    
    
    grDevices::dev.off()
    write.table(network_df, file=paste0(file.prefix, ".txt"), row.names = F, sep="\t")
    
  }
  
  
  
  return(list(igraph=path_mod, network_data=network_df, path_comm_genes=path_comm_genes, sim_coeff=set_sim_df))
  
}
