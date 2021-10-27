#' Caluclate cross-talks between pathways
#' @param pathway_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param wight an ordered weight vertex vector
#' @param membership an ordered vector of vertex membership to a topological community
#' @param mc_cores_pct numebr of threads to be used to calculate cross talk between different communities
#' @param mc_cores_perm numebr of threads to be used to calcualte permutation
#' @param mc_cores_tm number of threads to be used to calcualte TM-PCT on different communities combination
#' @param k number of permutation of the adjacency matrix
#' @import parallel
#' @import gtools
#' @importFrom kit funique
#' @export

TM_PCT <- function (pathway_list, gene_network_adj, membership = attr_v$comm_id,
                    mc_cores_pct = 2, mc_cores_perm = 2, mc_cores_tm = 2,
                    weight = attr_v$Sp,
                    k=9) {

  gene_network_adj <- sign(gene_network_adj)
  names(weight) <- rownames(gene_network_adj)
  perm_list <- mclapply(1:k, function(x) {
    tmp <- matrix(as.numeric(gene_network_adj), ncol = ncol(gene_network_adj),
                  dimnames = list(sample(rownames(gene_network_adj), nrow(gene_network_adj))))
    colnames(tmp) <- rownames(tmp)
    return(tmp)
  }, mc.cores = mc_cores_perm)


  unqMem <- unique(membership)
  mem_ptwL <- mclapply(unqMem, function(k) {
    tmp <- lapply(pathway_list, function(x) x <- x[x %in% rownames(gene_network_adj)[membership == k]])
    tmp <- tmp[which(lengths(tmp) != 0)]
  }, mc.cores = mc_cores_tm
  )

  names(mem_ptwL) <- unqMem
  mem_ptwL <- mem_ptwL[which(lengths(mem_ptwL) != 0)]

  comb_mem <- permutations(n=length(names(mem_ptwL)), r=2, v = names(mem_ptwL), repeats.allowed = F)
  comb_mem <- kit::funique(t(apply(comb_mem, 1, sort)))

  xx <- mclapply(1:nrow(comb_mem),  function (x) {

    idx_1 <- as.character(rownames(gene_network_adj)[membership == comb_mem[x,1]])
    idx_2 <- as.character(rownames(gene_network_adj)[membership == comb_mem[x,2]])

    tmp <- gene_network_adj[idx_1, idx_2, drop = F]

    return( tmp)
  }, mc.cores = mc_cores_pct
  )
  idx <- which(mclapply(xx, function(n) sum(n, na.rm = T), mc.cores = mc_cores_pct) !=0)
  if(length(idx) == 0) {
    print("No cambination of TM available to calculate PCT")
    return("No cambination of TM available to calculate PCT")
  } else {
    xx2 <- xx[idx]
    comb_mem <- comb_mem[idx,]


    res <- parallel::mclapply(1:nrow(comb_mem), function(j) {
      ptw_listCC1 <- mem_ptwL[[comb_mem[j,1]]]
      ptw_listCC2 <- mem_ptwL[[comb_mem[j,2]]]

      comb_p <- expand.grid(ptwCC1 = names(ptw_listCC1),
                            ptwCC2 = names(ptw_listCC2), stringsAsFactors = F)
      comb_p <- comb_p[which(comb_p$ptwCC1 != comb_p$ptwCC2),]

      xxCC <- mclapply(1:nrow(comb_p),  function (x) {
        idx_1 <- as.character(ptw_listCC1[[comb_p[x,1]]])
        idx_2 <- as.character(ptw_listCC2[[comb_p[x,2]]])

        tmp <- gene_network_adj[idx_1, idx_2, drop = F]

        return( tmp)
      }, mc.cores = mc_cores_pct)


      idxCC <- which(mclapply(xxCC, function(n) sum(n, na.rm = T), mc.cores = mc_cores_pct) !=0)
      if (length(idxCC) != 0) {
        xx2CC <- xxCC[idxCC]
        comb_p <- comb_p[idxCC,]
        xxPL <- mclapply(1:k, function(n) mclapply(1:nrow(comb_p),  function (x) {
          idx_1 <- as.character(ptw_listCC1[[comb_p[x,1]]])
          idx_2 <- as.character(ptw_listCC2[[comb_p[x,2]]])

          tmp <- perm_list[[n]][idx_1, idx_2, drop = F]

          return( tmp)
        }, mc.cores = mc_cores_pct),  mc.cores = mc_cores_perm)

        pct <- mclapply(1:length(xx2CC), function(m) {

          gene_network_adj_ijW <- t(weight[rownames(xx2CC[[m]])]) %*% as.matrix(xx2CC[[m]])
          gene_network_adj_ijW <- gene_network_adj_ijW %*% weight[colnames(xx2CC[[m]])]
          gene_network_adj_ijW <- data.frame(commID_1 = comb_mem[j,1],
                                             pathway_1 = comb_p[m,1],
                                             commID_2 = comb_mem[j,2],
                                             pathway_2 = comb_p[m,2],
                                             pct = gene_network_adj_ijW,
                                             ngenes_pathway1 = length(ptw_listCC1[[comb_p[m,1]]]),
                                             ngenes_pathway2 = length(ptw_listCC2[[comb_p[m,2]]]),
                                             nlink = sum(xx2CC[[m]]),
                                             stringsAsFactors = F)
        },   mc.cores = mc_cores_pct)
        pct <- do.call(rbind, pct)

        p_list <- mclapply(1:k, function(j) {
          tmp <-  mclapply(1:length(xxPL[[j]]), function(m) {

            gene_network_adj_ijW <- t(weight[rownames(xxPL[[j]][[m]])]) %*% as.matrix(xxPL[[j]][[m]])
            gene_network_adj_ijW <- gene_network_adj_ijW %*% weight[colnames(xxPL[[j]][[m]])]
            gene_network_adj_ijW <- data.frame(pathway_1 = comb_p[m,1],
                                               pathway_2 = comb_p[m,2],
                                               pct = gene_network_adj_ijW,
                                               stringsAsFactors = F)
          },   mc.cores = mc_cores_pct)
          ans <- do.call(rbind, tmp)

          return(ans)
        },  mc.cores = mc_cores_perm
        )
        p_list <- mclapply(p_list, function(x) {
          rownames(x) <- paste(x$pathway_1, x$pathway_2, sep = "|")
          x <- x[, "pct", drop = F]
          return(x)
        },  mc.cores = mc_cores_perm)
        tmp <- pct[,"pct", drop=F]
        rownames(tmp) <- paste(pct$pathway_1, pct$pathway_2, sep = "|")
        p_list <- c(list(tmp), p_list)

        p_val <- calc_p(p_list)
        out <- pct
        out$p_value <- p_val

        return(list(res = out,
                    p_list = p_list))
      } else {
        ans <- NULL
        return(ans)
      }





    }, mc.cores = mc_cores_tm)

    res <- res[which(lengths(res)!=0)]
    res1 <- mclapply(1:length(res), function(k) {
      tmp <- res[[k]]$res
      return(tmp)
    }, mc.cores = mc_cores_cc)
    res1 <- do.call(rbind, res1)
    perm <- mclapply(1:length(res), function(k) {
      tmp <- res[[k]]$p_list
      return(tmp)
    }, mc.cores = mc_cores_cc)

    res1$eFDR <- eFDR(real_values = as.vector(unlist(res1$pct)),
                      all_values = as.vector(unlist(perm)))

    res <- list(comm_pathway_list = mem_ptwL, TM_PCT_res = res1)
    return(res)
  }

}
