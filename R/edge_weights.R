#' Function for the assignment of edges weights in order to better visualize the communities
#' @param community list with the mandatory field "membership"
#' @param network pathway network
#' @param weight.within value to weight the attraction between two vertices of the same community
#' @param weight.between value to weight the attraction between two vertices of two distinct communities
#' @return vector of edge weights
#' @import igraph
#' @export
#'
edge_weights <- function(community, network, weight.within = 100, weight.between = 1) {
  bridges <- igraph::crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights)
}
