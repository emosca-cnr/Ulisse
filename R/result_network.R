#' Function to compare networks of PCT or CCC using Jaccard index on the edges and/or Structural similarity
#' @description This functionis used to calculate Jaccard index and/or Structural similarity to compare multiple PCTN or CCCN
#' @details This function takes as an input a list of networks. Then it calcuates Jaccard index on the edges and/or 
#' Structural similarityes between all pairs
#' @param res_list a named list of networks
#' @param type which measure to calcualte, can be either "jaccard", "strcutural" o both
#' @return The function returns a matrix for each mesure calcualted. The matrix has rows and column numbers equal to the
#' number of elements in `res_list` and named after the names of the elements in it
#' @import igraph
#' @export

result_network <- function(res_list, type = c("jaccard", "structural")) {
  
  to_use <- lapply(res_list, function(x) x <- V(x)$name)
  inters.to_use <- Reduce(union, to_use)
  
  res_list <- lapply(res_list, function(x) {
    to_add <- inters.to_use[!inters.to_use %in% V(x)$name]
    x <- add_vertices(x, nv = length(to_add), name = to_add)
    return(x)
  })
  
  if("jaccard" %in% type) {
    adj.res_list <- lapply(res_list, function(x) {
      x <- as_adjacency_matrix(x)
      x <- x[inters.to_use, inters.to_use]
      return(x)
    })
    
    res_out_j <- data.frame(matrix(data = 0, nrow = length(res_list), ncol = length(res_list),
                                 dimnames = list(names(res_list), names(res_list))), stringsAsFactors = F)
    for (i in 1:length(res_list)) {
      for (j in 1:length(res_list)) {
        if(i == j) {
          res_out[i, j] <- 1
        } else if (i > j) {
          tmp.1 <- adj.res_list[[i]]
          tmp.2 <- adj.res_list[[j]]
          tot <- tmp.1 + tmp.2
          num <- sum(tot == 2)
          den <- sum(tot > 0)
          jac <- num/den
          
          res_out_j[i, j] <- jac
          res_out_j[j, i] <- jac
          
        }
      }
      
    }
  }
  
  if("structural" %in% type ) {
    
    res_out_s <- data.frame(matrix(data = 0, nrow = length(res_list), ncol = length(res_list),
                                   dimnames = list(names(res_list), names(res_list))), stringsAsFactors = F)
    for (i in 1:length(res_list)) {
      for (j in 1:length(res_list)) {
        if(i == j) {
          res_out[i, j] <- 1
        } else if (i > j) {
          tmp.1 <- res_list[[i]]
          tmp.2 <- res_list[[j]]
          sim <- D(g = tmp.1, h = tmp.2, w1 = 0.45, w2 = 0.45, w3 = 0.1)
          
          res_out_s[i, j] <- sim
          res_out_s[j, i] <- sim
          
        }
      }
    }
    
  }
  
  if("jaccard" %in% type) {
    if("structural" %in% type) {
      return(list(Jaccard = res_out_j,
                  Structural = res_out_s))
      } else {
        return(res_out_j)
      }
    
    } else {
      return(res_out_s)
    }
}
