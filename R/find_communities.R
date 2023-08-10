#' Find topological communities
#' @param g graph
#' @param e.weights edge weights
#' @param v.weights vertex weights
#' @param verbose TRUE/FALSE
#' @param methods one or more of "fastgreedy", "labprop", "walktrap", "eigen", "multilev", "infomap"
#' @return list od community objects
#' @import igraph
#' @export
find_communities <- function(g, e.weights=NULL, v.weights=NULL, verbose=TRUE, methods=c("fastgreedy", "labprop", "walktrap", "eigen", "multilev", "infomap")){


  comm_list <- vector('list', 6)
  names(comm_list) <- methods
  comm_info <- data.frame(algorithm=methods,
                          modularity=NA,
                          n=NA,
                          stringsAsFactors = FALSE
  )

  if("fastgreedy" %in% methods){
    if(verbose)
      print("fastgreedy")
    idx <- which(grepl("fastgreedy", methods))
    comm_list[[idx]] <- igraph::fastgreedy.community(g, weights = e.weights)
    comm_info$modularity[idx] <- igraph::modularity(g, comm_list[[idx]]$membership, weights = e.weights)
    comm_info$n[idx] <- max(comm_list[[idx]]$membership)

  }
  if("labprop" %in% methods){
    if(verbose)
      print("labprop")
    idx <- which(grepl("labprop", methods))
    comm_list[[idx]] <- igraph::label.propagation.community(g, weights = e.weights)
    comm_info$modularity[idx] <- igraph::modularity(g, comm_list[[idx]]$membership, weights = e.weights)
    comm_info$n[idx] <- max(comm_list[[idx]]$membership)
  }
  if("walktrap" %in% methods){
    if(verbose)
      print("walktrap")
    idx <- which(grepl("walktrap", methods))
    comm_list[[idx]] <- igraph::walktrap.community(g, weights = e.weights)
    comm_info$modularity[idx] <- igraph::modularity(g, comm_list[[idx]]$membership, weights = e.weights)
    comm_info$n[idx] <- max(comm_list[[idx]]$membership)
  }

  if("eigen" %in% methods){
    if(verbose)
      print("eigen")
    idx <- which(grepl("eigen", methods))
    comm_list[[idx]] <- igraph::leading.eigenvector.community(g, weights = e.weights)
    comm_info$modularity[idx] <- igraph::modularity(g, comm_list[[idx]]$membership, weights = e.weights)
    comm_info$n[idx] <- max(comm_list[[idx]]$membership)
  }

  if("multilev" %in% methods){
    if(verbose)
      print('multilev')
    idx <- which(grepl("multilev", methods))
    comm_list[[idx]] <- igraph::multilevel.community(g, weights = e.weights)
    comm_info$modularity[idx] <- igraph::modularity(g, comm_list[[idx]]$membership, weights = e.weights)
    comm_info$n[idx] <- max(comm_list[[idx]]$membership)
  }

  if("infomap" %in% methods){
    if(verbose)
      print('infomap')
    idx <- which(grepl("infomap", methods))
    comm_list[[idx]]  <- igraph::infomap.community(g, e.weights = e.weights, v.weights = v.weights)
    comm_info$modularity[idx] <- igraph::modularity(g, comm_list[[idx]]$membership, weights = e.weights)
    comm_info$n[idx] <- max(comm_list[[idx]]$membership)
  }


  return(list(comm=comm_list, info=comm_info))

}
