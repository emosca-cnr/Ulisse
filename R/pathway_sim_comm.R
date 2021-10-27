#' Pathway similarity and communities
#' @param x named vector of pathway scores
#' @param gs_list gene set list
#' @param method 'overlap' or 'jaccard'
#' @param coeff threshold for the similarity score between two gene sets
#' @param all_gs TRUE/FALSE
#' @param comm_method one among 'auto', 'fastgreedy' and 'multilev'
#' @param set_sim_df optional data frame with three columns 'set1', 'set2' and 'sim'
#' @return list of two data frames containing pathway network and vertex similarity
#' @import igraph
#' @importFrom  RColorBrewer brewer.pal
#' @export
#'
pathway_sim_comm <- function(x, gs_list, method=c('overlap', 'jaccard'), coeff=c(0.5, 0.1), all_gs=TRUE, comm_method=c('auto', 'fastgreedy', 'multilev'), set_sim_df=NULL){

  method <- match.arg(method)
  comm_method <-  match.arg(comm_method)
  if(method == "overlap"){
    coeff <- coeff[1]
  }
  if(method == 'jaccard'){
    corff <- coeff[2]
  }

  gs_list_subset <- gs_list[names(gs_list) %in% names(x)]


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
  set_sim_df <- set_sim_df[set_sim_df$sim >= coeff, ]
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


  out <- data.frame(id=names(x), score=x, comm=0, size=1, stringsAsFactors = FALSE)
  out$comm[out$id %in% path_mod_comm$names]<- path_mod_comm$membership[match(out$id[out$id %in% path_mod_comm$names], path_mod_comm$names)]
  out$size[out$id %in% path_mod_comm$names] <- table(path_mod_comm$membership)[out$comm[out$id %in% path_mod_comm$names]]

  return(list(comm=out, sim_coeff=set_sim_df))

}
