#' Caluclate cross-talks between pathways
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param wight an ordered weight vertex vector
#' @param membership an ordered vector of vertex membership to a topological community
#' @param cluster numebr of threads to be used
#' @import parallel
#' @import doParallel
#' @export


TM_PCT <- function (pathway_list, gene_network_adj, membership = attr_v$comm_id, 
                    weight = attr_v$Sp, cluster = 4) {
  comb_mem <- expand.grid(comm_1 = unique(membership),
                          comm_2 = unique(membership), stringsAsFactors = F)
  comb_mem <- data.frame(comb_mem[ !duplicated(apply(comb_mem, 1, sort), MARGIN = 2), ], 
                         stringsAsFactors = F)
  comb_mem <- comb_mem[-which(comb_mem$comm_1 == comb_mem$comm_2),]
  
  gene_network_adj <- data.frame(sign(gene_network_adj), stringsAsFactors = F)
  colnames(gene_network_adj) <- rownames(gene_network_adj)
  
  res <- parallel::mclapply(1:nrow(comb_mem), function(x) {
    
    idxCC1 <- which(membership == comb_mem[x,1])
    # gene_network_adjCC1 <- gene_network_adj[idxCC1, idxCC1]
    ptw_listCC1 <- lapply(pathway_list, function(x) x <- x[x %in% colnames(gene_network_adj)[idxCC1]])
    ptw_listCC1 <- ptw_listCC1[which(lengths(ptw_listCC1) != 0)]
    weightCC1 <- weight[idxCC1]
    
    idxCC2 <- which(membership == comb_mem[x,2])
    #gene_network_adjCC2 <- gene_network_adj[idxCC2, idxCC2]
    ptw_listCC2 <- lapply(pathway_list, function(x) x <- x[x %in% colnames(gene_network_adj)[idxCC2]])
    ptw_listCC2 <- ptw_listCC2[which(lengths(ptw_listCC2) != 0)]
    weightCC2 <- weight[idxCC2]
    
    comb_p <- expand.grid(ptwCC1 = names(ptw_listCC1),
                          ptwCC2 = names(ptw_listCC2), stringsAsFactors = F)
    comb_p <- data.frame(comb_p[ !duplicated(apply(comb_p, 1, sort), MARGIN = 2), ], 
                         stringsAsFactors = F)
    if(length(which(comb_p$ptwCC1 == comb_p$ptwCC2))>0) {
      comb_p <- comb_p[-which(comb_p$ptwCC1 == comb_p$ptwCC2),]
    }
    
    #n <- paste(comb_p$ptw1, comb_p$ptw2, sep = "|")
    if (nrow(comb_p)>0) {
      xx <- apply(comb_p, 1, function (j) {
        idx_1 <-  which(rownames(gene_network_adj) %in% ptw_listCC1[[j[1]]])
        idx_2 <-  which(rownames(gene_network_adj) %in% ptw_listCC2[[j[2]]])
        return( gene_network_adj[idx_1, idx_2])
      }
      )
      names(xx) <- paste(comb_p[,1], comb_p[,2], sep = "|")
      idx <- which(lapply(xx, function(n) sum(n)) !=0)
      xx2 <- xx[idx]
      
      
      if(length(xx2) > 0) {
        pct <- mclapply(names(xx2), function(y) {
          p_i <- unlist(lapply(strsplit(y, "|", fixed = T), "[[", 1))
          p_j <- unlist(lapply(strsplit(y, "|", fixed = T), "[[", 2))
          
          idx_i <- which(rownames(gene_network_adj) %in% ptw_listCC1[[p_i]])
          idx_j <- which(rownames(gene_network_adj) %in% ptw_listCC2[[p_j]])
          i_genes <- rownames(gene_network_adj)[idx_i]
          j_genes <- rownames(gene_network_adj)[idx_j]
          gene_network_adj_ij <- as.matrix(xx2[[y]]) 
          if (length(which(i_genes %in% j_genes))>1) {
            
            s <- i_genes[which(i_genes %in% j_genes)]
            sidx_i <- which(i_genes %in% s)
            sidx_j <- which(j_genes %in% s)
            
            gene_network_adj_ij <- gene_network_adj_ij[-sidx_i, -sidx_j]
            i_genes <- i_genes[-sidx_i]
            j_genes <- j_genes[-sidx_j]
            idx_i <- idx_i[-sidx_i]
            idx_j <- idx_j[-sidx_j]
            
          }
          #print(y)
          
          gene_network_adj_ijW <- t(weight[idx_i]) %*% gene_network_adj_ij
          gene_network_adj_ijW <- gene_network_adj_ijW %*% weight[idx_j]
          gene_network_adj_ijW <- data.frame(pathway_1 = p_i,
                                             pathway_2 = p_j,
                                             pct = gene_network_adj_ijW,
                                             stringsAsFactors = F)
          colnames(gene_network_adj_ijW)[1:2] <- paste0(rep("comm_", 2), comb_mem[x,])
          return(gene_network_adj_ijW)
        },   mc.cores = cluster)
        ans <- do.call(rbind, pct)
        ans <- list(pathway_list1=ptw_listCC1, patwhay_list2 = ptw_listCC2,
                    TM_PCT = ans)
        names(ans)[1:2] <- paste0(rep("comm_", 2), comb_mem[x,], "_pathway_list")
        return(ans)
      } 
    }
    
    
  }, mc.cores = cluster)
  
  names(res) <- paste(comb_mem[,1], "VS", comb_mem[,2], sep = "_")
  return(res) 
}
