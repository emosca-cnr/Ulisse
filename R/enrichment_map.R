#' Enrichment map
#' @description This functions calculates similarities between gene-sets and plot a resulting enrichment map
#' @details enrichment_map() function calculates similarities between each gene-set pair by using `method` metric. Subsequently,
#'  these similarities are filtered to maintain the ones >= `coeff`. `comm_method` algorithm is then used to identify communities, 
#'  which may be filtered to plot only the ones composed by at least `min_comm_size` gene-sets. Enrichment map is then plotted by using
#'  `ggraph` package.
#' @param gs_score named vector of pathway scores
#' @param gs_list gene set list
#' @param method 'overlap' or 'jaccard'
#' @param coeff threshold for the similarity score between two gene sets
#' @param all_gs TRUE/FALSE, indicating if all gene-sets should be considered or only the ones with similarity score >= `coeff`
#' @param comm_method community algorithm to be used; available algorithms are c("fastgreedy", "labprop", "walktrap", "eigen", "multilev", "infomap")
#'  see `Ulisse::find_communities()` for further details
#' @param min_comm_size minimum size of a community to be considered in the enrichment map. If `min_comm_size = 1` and `all_gs = TRUE`,
#'  then all gene-sets are displayed
#' @param gs_list_size named vector with size f the gene-sets in `gs_list`
#' @param set_sim_df optional, data.frame with three columns 'set1', 'set2' and 'sim'
#' @param file_prefix prefix used to save similarity data.frame and enrichment map plot (if `save_plot = TRUE`)
#' @param layout optional layout matrix, composed of two columns. If not provided, the layout is calculated by using
#'  `igraph::layout_with_fr()` function, together with `weight_within` and `weight_between` (via `Ulisse::edge_weights()` function)
#' @param weight_within value to weight the attraction between two vertices of the same community
#' @param weight_between value to weight the attraction between two vertices of two distinct communities
#' @param pal_community palette to be used to color communities, should be a named vector with a color for each community. If not
#'  provided, `pals::alphabet2()` is used instead
#' @param pal_score palette to be used to color vertices according to their `score`. If not provided, `pals::brewer.greens(n = 3)` is used
#' @param top_ptw variable used to select top gene-sets to display as a description of the communities. 
#'  Should be either `"score"` (to display top score gene-set for each community) or `"n_genes"` (to display gene-set with higher number of genes)
#' @param wrap value used to wrap long gene-sets names as a description of communities. See `stringr::str_wrap()` for further details
#' @param label_fontsize font size of the label and description of the communities. If two values are provided,
#'  the fist is used for label, the second for description. See `ggforce::geom_mark_rect()` for further details.
#' @param e_color edge color
#' @param e_alpha edge transparency
#' @param e_range range of widths values used for plotting edges proportional to similarity score
#' @param v_stroke width of the stoke of the vertices, which is colored according to communities
#' @param v_range range of dimensions of the vertices, which is proportional to the number of genes
#' @param save_plot TRUE/FALSE, if the plot should be saved by using `file_prefix` or returned
#' @param width jpeg image width
#' @param height jpeg image height
#' @param res jpeg image resolution
#' @param units jpeg image units for width and height
#' @return The function returns a list containing: 
#' \itemize{
#'  \item igraph = network object used for plotting the enrichment map
#'  \item network_data = data.frame with layout (X1,X2 columns), name of the gene-sets together with community ids ("comm_id"), 
#'   score and number of genes ("n_genes")
#'  \item path_comm_genes = list composed by the genes present in each community
#'  \item sim_coeff = data.frame with the similarity score calculated between each gene-set pair
#'  \item plot = enrichment map plot obtained by using `ggraph` package functions. Only if `save_plot = FALSE`
#' }
#' @import igraph pals stats graphics
#' @importFrom utils write.table
#' @importFrom grDevices dev.off jpeg
#' @import ggplot2
#' @import ggraph
#' @import pals
#' @importFrom ggforce geom_mark_rect
#' @importFrom ggnewscale new_scale_fill
#' @import stringr
#' @export
#'
enrichment_map <- function(gs_score, gs_list, method=c('overlap', 'jaccard'),
                           coeff=NULL, all_gs=TRUE, comm_method='fastgreedy', min_comm_size=2,
                           gs_list_size=NULL, set_sim_df=NULL, file_prefix='en_map', 
                           layout = NULL, weight_within = 4, weight_between = 1,
                           pal_community = NULL, pal_score = NULL, top_ptw = "score",
                           wrap = 15, label_fontsize = c(6, 5), 
                           e_color = "gray65", e_alpha = 0.5, e_range = c(1,4), 
                           v_stroke = 0.5, v_range = c(2,8), 
                           save_plot=TRUE, width=200, 
                           height=200, res=300, units = "mm"){
  #method=c('overlap')
  #similarity method and coeff----------------
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
    gs_list <- gs_list[names(gs_list) %in% names(gs_score)]
    gs_list_size <-lengths(gs_list)
  }
  
  #calculating similarities-----------------------
  if(is.null(set_sim_df)){
    cat("Calculating gene set similarities... ")
    set_sim_df <- data.frame(t(utils::combn(names(gs_score), 2)), stringsAsFactors = FALSE) #build table of pair occurrences
    colnames(set_sim_df)[1:2] <- c('set1', 'set2')
    cat( nrow(set_sim_df), 'similarities found\n')
    set_sim_df$sim <- sapply(1:nrow(set_sim_df), function(x) {
      s.1 <- as.character(set_sim_df[x, 1])
      s.2 <- as.character(set_sim_df[x, 2])
      out <- calc_set_similarity(gs_list[[s.1]], 
                                 gs_list[[s.2]], method=method)
      return(out)
    }, simplify = T)
    
    write.table(set_sim_df, paste0(file_prefix, "_sim_table.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
    
    cat("done.\n")
  }else{
    cat("Using given gene set similarities\n")
    cat("Keep only gene sets for which similarity values are avilable\n")
    gs <- unique(c(set_sim_df[, 1], set_sim_df[, 2]))
    gs_list <- gs_list[names(gs_list) %in% gs]
    gs_score <- gs_score[names(gs_score) %in% names(gs_list)]
  }
  
  
  
  cat("Summary of gene set similarities (", nrow(set_sim_df),"):\n")
  print(summary(set_sim_df$sim))
  
  #building pathway network + communities------------------------------
  path_mod <- igraph::graph.data.frame(set_sim_df[set_sim_df$sim >= coeff, ], 
                                                        directed=FALSE)#build network
  if(length(igraph::V(path_mod)) < length(gs_score) & all_gs){
    to_add <- names(gs_score)
    to_add <- to_add[!to_add %in% V(path_mod)$name]
    l <- length(to_add)
    path_mod <- igraph::add_vertices(path_mod, l, 
                                     attr = list(name=to_add)) #reintroduce vertexes excluded
  } #reintroduce vertex at degree 0
  cat("Enrichment Map: V=", length(V(path_mod)), ", E=", length(E(path_mod)), "\n")
  
  cat("Community detection...\n")
  path_mod_comm <- find_communities(path_mod, e.weights = igraph::E(path_mod)$sim, methods = comm_method)
  print(path_mod_comm$info)
  
  path_mod_comm <- path_mod_comm$comm[[comm_method]]
  
  V(path_mod)$comm_id <- as.character(path_mod_comm$membership)
   
  if(min_comm_size > 1){ #filtering network by community size. Removing communities smaller that `min_comm_size`
    cat("Excluding communities smaller than ", min_comm_size, "\n")
    keep_comm <- sizes(path_mod_comm)
    keep_comm <- names(keep_comm)[keep_comm >= min_comm_size]
    cat("Kept:", length(keep_comm), "\n")
    path_mod <- igraph::induced_subgraph(path_mod, igraph::V(path_mod)$name[igraph::V(path_mod)$comm_id %in% keep_comm])
    cat("Enrichment Map: V=", length(V(path_mod)), ", E=", length(E(path_mod)), "\n")
  }
  
  V(path_mod)$score <- gs_score[V(path_mod)$name] #pathway score based on ORA res
  V(path_mod)$size <- gs_list_size[V(path_mod)$name]
  ### genes of each pathway community
  path_comm_genes <- split(V(path_mod)$name, V(path_mod)$comm_id)
  path_comm_genes <- lapply(path_comm_genes, function(pathway_vector) sort(unique(unlist(gs_list[names(gs_list) %in% pathway_vector]))))
  
  #layout---------------------------------
  if(is.null(layout)){
    layout <- igraph::layout_with_fr(path_mod, weights = edge_weights(list(membership=V(path_mod)$comm_id), path_mod, 
                                                                      weight.within = weight_within, weight.between = weight_between))
  }
  
  #preparing tables for returning-----------------
  network_df <- data.frame(layout, name=V(path_mod)$name, comm_id=V(path_mod)$comm_id, 
                           score=V(path_mod)$score, stringsAsFactors = F) ##all
  network_df$n_genes <- gs_list_size[network_df$name]
  tmp <- split(network_df, network_df$comm_id)
  n <- names(tmp)
  tmp <- unlist(lapply(tmp, function(tab) {
    tab <- tab[order(tab[, top_ptw], decreasing = F),]
    tab <- tab$name[1]
    tab <- gsub("_", " ",  tab)
    return(tab)
  }))
  
  names(tmp) <- n
  V(path_mod)$comm_id2 <- tmp[as.character(V(path_mod)$comm_id)]
  
  #palettes---------------------------
  if(is.null(pal_community)) {
    pal_community <- pals::alphabet2(n = length(n))
    names(pal_community) <- n
  } else {
    pal_community <- pal_community[n]
  }
  
  
  if(is.null(pal_score)) {
    pal_score <- pals::brewer.greens(n = 3)
  } 
  
  
  #plot--------------------------------
  
  #set pals + alpha values
  #names of the pathway with most genes/highest score
  #community labels
  p <- ggraph(path_mod, layout) +
    theme_graph() +
    ggforce::geom_mark_rect(aes(x = x, y = y, label = comm_id, 
                                description = stringr::str_wrap(comm_id2, wrap), 
                                fill = comm_id, color = comm_id), 
                            show.legend = FALSE, label.fontsize = label_fontsize, con.cap = unit(1, "mm"),
                            label.margin = margin(1, 0, 1, 0, "mm")) +
    scale_fill_manual(values = pal_community) +
    scale_color_manual("Community", values = pal_community) +
    geom_edge_link(aes(edge_width = sim), color = e_color, alpha = e_alpha) +
    scale_edge_width("Similarity", range = e_range) +
    ggnewscale::new_scale_fill() +
    geom_node_point(aes(fill = score, size = size, color = comm_id), shape = 21, stroke = v_stroke) +
    scale_fill_gradientn(colours = pal_score) +
    scale_size("# genes", range = v_range) +
    guides(color = guide_legend(override.aes = list(shape = 19, size = 3))) #+
    # theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
    #       legend.key.size = unit(5, 'mm'))
  
  if(save_plot) {
    jpeg(paste0(file_prefix, ".jpeg"), width = width, height = height, 
         res = res, units = units)
  }
  print(p)
  if(save_plot) {
    dev.off()
  }
  
  if(save_plot) {
    return(list(igraph=path_mod, network_data=network_df, path_comm_genes=path_comm_genes, sim_coeff=set_sim_df))
  } else {
    return(list(igraph=path_mod, network_data=network_df, path_comm_genes=path_comm_genes, sim_coeff=set_sim_df, plot = p))
  }
  
  
}

