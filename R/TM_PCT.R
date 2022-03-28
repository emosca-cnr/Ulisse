#' Function to caluclates cross-talks between pathways of differente gene communities
#' @description Calculates cross-talks between pathways composed by gene of different gene communities
#' @details The function uses the `membership` data to reogranize the `pathway_list` as to obtain a list of 
#' pathway list for each gene communities. In detail, gene will be grouped into communities and then grouped into pathways.
#' At this point, the function calculates CT between each pathway pairs of different 
#' communities that shows at least a link.
#' This approach is applied also to the permuted version of the `gene_network_adj`. 
#' CT and permutd CT are then used to calcualte `p-value` and `FDR`
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param weight weights of the genes in gene_network_adj
#' @param membership membership of the gene in gene_network_adj to topological communities
#' @param mc_cores_pct numebr of threads to be used to calculate cross talk 
#' @param mc_cores_perm numebr of threads to be used to calcualte permutation
#' @param mc_cores_tm number of threads to be used to calcualte TM-PCT on different communities combination
#' @param k number of permutation of the adjacency matrix
#' @return The function returns a list of two object:
#' \enumerate{
#'  \item comm_pathway_list: the community pathway list used for TM-PCT calculation. Each list is named with the name 
#'  of the gene community
#'  \item {TM_PCT_res}: a data frame with
#' \itemize{
#'  \item commID_1, commID_2: community to which belongs the genes in `pathway_1` and `pathway_2`, respectively
#'  \item pathway_1, pahtway_2: the pathway pair considered
#'  \item score: TM-PCT score
#'  \item ngenes_pathway_1, ngenes_pathway_2: number of genes involved in `pathway_1` and `pathway_2`, respectively
#'  \item n_link: number of links between the pathways considered
#'  \item gene_pathway1, gene_pathway2: gene involevd in the CT in `pathway_1` and `pathway_2`, respectively
#'  \item p_value: empirical p-value calculated by using the permutation approcah
#'  \item eFDR: empirical FDR
#'  }
#'  }
#' @examples  
#'  ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#'  adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#'  ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#'  wgt <- rep(1, 6)
#'  memb <- c(1, 1, 2, 2, 3, 3)
#'  pct <- TM_PCT(pathway_list = ptw_list, gene_network_adj = adj, weight = wgt, membership = memb, 
#'                 mc_cores_tm = 1, mc_cores_pct = 1, mc_cores_perm = 1, k = 9)
#' @import parallel
#' @import gtools
#' @import kit
#' @import stringi
#' @export

TM_PCT <- function (pathway_list, gene_network_adj, membership, 
                    mc_cores_pct = 2, mc_cores_perm = 2, mc_cores_tm = 2,
                    weight, 
                    k=9
) {
  
  gene_network_adj <- sign(gene_network_adj)
  names(weight) <- rownames(gene_network_adj)
  names(membership) <- rownames(gene_network_adj)
  weight <- weight[weight !=0]
  gene_network_adj <- gene_network_adj[names(weight), names(weight)]
  membership <- membership[names(weight)]
  
  unqMem <- unique(membership)
  mem_ptwL <- mclapply(unqMem, function(z) {
    tmp <- lapply(pathway_list, function(x) x <- x[x %in% rownames(gene_network_adj)[membership == z]])
    tmp <- tmp[lengths(tmp) != 0]
  }, mc.cores = mc_cores_tm
  )
  
  names(mem_ptwL) <- unqMem 
  mem_ptwL <- mem_ptwL[lengths(mem_ptwL) != 0]
  
  comb_mem <- permutations(n=length(names(mem_ptwL)), r=2, v = names(mem_ptwL), repeats.allowed = F)
  comb_mem <- funique(t(apply(comb_mem, 1, sort)))
  
  xx <- mclapply(1:nrow(comb_mem),  function (x) {
    
    g.1 <- names(membership)[membership %in% comb_mem[x,1]]
    g.2 <- names(membership)[membership %in% comb_mem[x,2]]
    idx_1 <- as.character(g.1[!g.1 %in% g.2])
    idx_2 <- as.character(g.2[!g.2 %in% g.1])
    
    tmp <- gene_network_adj[idx_1, idx_2, drop = F]
    
    return( tmp) 
  }, mc.cores = mc_cores_pct
  )
  idx <- which(mclapply(xx, function(n) sum(n, na.rm = T), mc.cores = mc_cores_pct) !=0)
  if(length(idx) == 0) {
    print("No cambination of TM available to calculate PCT") 
    return("No cambination of TM available to calculate PCT")
  } else {
    xx <- xx[idx]
    comb_mem <- comb_mem[idx,,drop = F]
    
    perm_list <- mclapply(1:k, function(x) {
      tmp <- gene_network_adj
      rownames(tmp) <- sample(rownames(tmp), nrow(tmp))
      colnames(tmp) <- rownames(tmp)
      
      return(tmp)
    }, mc.cores = mc_cores_perm)
    
    
    res <- parallel::mclapply(1:nrow(comb_mem), function(j) {
      ptw_listCC1 <- mem_ptwL[[comb_mem[j,1]]]
      ptw_listCC2 <- mem_ptwL[[comb_mem[j,2]]]
      
      comb_p <- expand.grid(ptwCC1 = names(ptw_listCC1),
                            ptwCC2 = names(ptw_listCC2), stringsAsFactors = F)
      comb_p <- comb_p[comb_p$ptwCC1 != comb_p$ptwCC2,]
      
      xxCC <- mclapply(1:nrow(comb_p),  function (x) {
        idx_1 <- as.character(ptw_listCC1[[comb_p[x,1]]])
        idx_2 <- as.character(ptw_listCC2[[comb_p[x,2]]])
        tmp <- gene_network_adj[idx_1, idx_2, drop = F]
        
        return( tmp) 
      }, mc.cores = mc_cores_pct)
      
      
      idxCC <- which(mclapply(xxCC, function(n) sum(n, na.rm = T), mc.cores = mc_cores_pct) !=0)
      if (length(idxCC) != 0) {
        xxCC <- xxCC[idxCC]
        comb_p <- comb_p[idxCC,]
        
        xxPL <- mclapply(1:k, function(n) {
          tab <- perm_list[[n]]
          tmp <- mclapply(1:nrow(comb_p),  function (x) {
            idx_1 <- as.character(rownames(xxCC[[x]]))
            idx_2 <- as.character(colnames(xxCC[[x]]))
            tab.tmp <- tab[idx_1, idx_2, drop = F]
            return( tab.tmp) 
          }, mc.cores = mc_cores_pct)
        } ,  mc.cores = mc_cores_perm)
        
        allCC <- c(list(xxCC), xxPL)
        out <- mclapply(1:(k+1), function(n) {
          pct <- mclapply(1:length(xxCC), function(m) {
            tab <- allCC[[n]][[m]]
            g.1 <- rownames(tab)
            g.2 <- colnames(tab)
            
            gene_network_adj_ijW <- t(weight[g.1]) %*% tab
            gene_network_adj_ijW <- gene_network_adj_ijW %*% weight[g.2]
            
            row.col.idx <- which(tab == 1, arr.ind = T)
            row.n <- g.1[row.col.idx[,1]]
            row.n <- row.n[row.n %in% names(weight)]
            col.n <- g.2[row.col.idx[,2]]
            col.n <- col.n[col.n %in% names(weight)]
            
            tab.out <- array(data = NA, dim = 12)
            tab.out[1] <- comb_mem[j, 1]
            tab.out[2] <- comb_p[m,1]
            tab.out[3] <- comb_mem[j, 2]
            tab.out[4] <- comb_p[m,2]
            tab.out[5] <- gene_network_adj_ijW
            tab.out[6] <- length(row.n)
            tab.out[7] <- length(col.n)
            tab.out[8] <- sum(tab)
            tab.out[9] <- sum(weight[g.1])
            tab.out[10] <- sum(weight[g.2])
            tab.out[11] <- stri_c(row.n, collapse = ";")
            tab.out[12] <- stri_c(col.n, collapse = ";")
            
            
            return(tab.out)
          },   mc.cores = mc_cores_pct)
          
          pct <- do.call(rbind, pct)
          colnames(pct) <- c("commID_1", "pathway_1", "commID_2", "pathway_2", "pct", 
                             "ngenes_pathway1", "ngenes_pathway2", "nlink", "weight_pathway1", 
                             "weight_pathway2", "gene_pathway1", "gene_pathway2")
          return(pct)
        },  mc.cores = mc_cores_perm)
        
        
        return(out)
      } else {
        ans <- NULL
        return(ans)
      }
      
    }, mc.cores = mc_cores_tm)
    
    p_list <- mclapply(1:length(res), function(j) {
      tab_list <- res[[j]]
      out <- mclapply(tab_list, function(x) {
        tmp <- matrix(as.numeric(x[,"pct"]), ncol = 1)
        rownames(tmp) <- paste(x[,1], x[,2], sep = "|")
        return(tmp)
      }, mc.cores = mc_cores_perm)
      out.2 <- calc_p(out)
      return(out.2)
    }, mc.cores = mc_cores_pct)
    
    fdr_real <- lapply(res, function(x) {
      tmp <- x[[1]][, "pct"]
      return(tmp)
    })
    fdr_all <- unlist(lapply(res, function(x) {
      tmp <- lapply(x, function(j) j <- j[, "pct"])
    }))
    
    res_FDR <- lapply(1:length(fdr_real), function (x) {
      tmp <- eFDR(real_values = fdr_real[[x]], 
                  all_values = fdr_all)
      return(tmp)
    })
    out <- lapply(1:length(res), function(x) {
      tmp <- res[[x]][[1]]
      coln <- colnames(tmp)
      tmp <- cbind(tmp, p_list[[x]])
      tmp <- cbind(tmp, res_FDR[[x]])
      colnames(tmp) <- c(coln, "p_value", "FDR")
      return(tmp)
    })
    
    out <- do.call(rbind, out)
    out <- data.frame(out, stringsAsFactors = F)
    rownames(out) <- NULL
    
    res <- list(comm_pathway_list = mem_ptwL, TM_PCT_res = out)
    return(res)
  }
  
}
