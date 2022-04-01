#' Function to caluclate permutation of the number of links between two gene sets
#' @description `perm_link` is used to obtain permutatied number of links between two genesets of 
#' a defined number of gene
#' @details The function generates k permuted verion of two geneset of defined length by random sampling the 
#' gene network adjacency matrix. Then identify the subset of the adjacency matrix for each k geneset and 
#' calculates the number of links.
#' @param r length of the first geneset (rows)
#' @param c length of the second geneset (columns)
#' @param gene_network_adj adjacency matrix of the gene network used
#' @param core threads to be used for permutation 
#' @param k number of permutation
#' @param hash logical, indicating if random sampling should be done by using hashmap (see `sample.int()`)
#' @return The function a list k matrices, each one composed by the number of links in the permutated verion of the 
#' two genesets
#' @import parallel
#' @export

perm_link <- function(r, c, gene_network_adj, core, k, hash = T) {
  n.tot <- gene_network_adj@Dim[1]
  perm.1 <- mclapply(1:k, function(j) j <- sample.int(n.tot, r, useHash = hash),
                     mc.cores = core)
  perm.2 <- mclapply(1:k, function(j) j <- sample.int(n.tot, c, useHash = hash),
                     mc.cores = core)
  
  perm <- mclapply(1:k, function(j) {
    
    tmp.1 <- perm.1[[j]]
    tmp.2 <- perm.2[[j]]
    perm_tab <- gene_network_adj[tmp.1, tmp.2]
    perm_l <- matrix(sum(perm_tab))
    return(perm_l)
  },   mc.cores = core)
  return(perm)
  
}



