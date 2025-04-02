#' Integrated cellular cross talk analysis
#' @description The function is a wrapper of `single_integrated_cross_talk()` and calculates inter-intracellular communication on multiple clusters. 
#' @details The function takes as inputs takes as input the list of communication details of interest (`ccc_list`), 
#'  together with the cluster-list used for CCC calculation (`cl_list`), a list of gene sets lists (`gs_list`), one for each cluster/cell type,
#'  and the adjacency matrix of the biological network (`A`). 
#'  These inputs are formatted and passed to `single_integrated_cross_talk()` to calculate integrated cross-talk on each cluster pair provided. 
#' @param ccc_list communication details between cluster of interest
#' @param cl_list cluster list used for communication analysis (provided as `S_list` to `gs_cross_talk()`)
#' @param gs_list list of gene set lists, one for each cluster considered in the analysis
#' @param A adjacency matrix of the whole gene network considered (can be a sparseMatrix)
#' @param k number of permutations
#' @param bin_type can be either "number" (suggested), "interval" or "width". See `ggplot2::cut_interval()` for details.
#' @param cut_par number of bin to cut. If set to `NULL`, the function will search for the best cut in 2:15 
#' @param perm_link = c("degree", "simple") If the permutation of the link should be degree-conservative ("degree"), or random ("simple")
#' @param perm_weights = c("degree", "simple") If the permutation of the weights should be degree-conservative ("degree"), or random ("simple")
#' @param mc_cores_ct number of threads to be used for cross talk calculation
#' @param mc_cores_perm number of thread to be used in permutations
#' @return The function returns a data.frame with the result of inter-intracellular communication:
#' \itemize{
#'    \item ccc: cluster pair considered for integrated cross-talk analysis
#'    \item S1_name,S2_name: names of gene-set pairs considered
#'    \item c: cross-talk score
#'    \item S1_size,S2_size: number of genes in S1 and S2
#'    \item S1_s2_size,S2_s1_size: number of S1 genes interacting with S2, and vice versa
#'    \item dL: number of links between S1 and S2
#'    \item L: number of possible links between S1 and S2
#'    \item r_c: cross-talk saturation, calculated as `dL/L`
#'    \item u1,u2: sum of the gene weights in S1 and S2, respectively
#'    \itme S1,S2: list of interacting genes in S1 and S2, respectively
#'    \item s: cross-talk summary score
#'    \item pA,pU: p-values of the number of links (pA) and weights (pU)
#'    \item p: combined p-value
#'  }
#' @importFrom parallel mclapply
#' @importFrom igraph simplify graph_from_edgelist as_edgelist degree
#' @importFrom stringi stri_c
#' @importFrom reshape2 acast
#' @importFrom methods is as
#' @importFrom collapse ss get_elem rowbind
#' @export
integrated_cross_talk <- function(cl_list = cl_list, ccc_list = ccc_list, 
                          gs_list = ptw_list, A = adj.ppi, k = 9, 
                          perm_link = "degree", perm_weights = "degree", cut_par = 3,
                          mc_cores_perm = 1, mc_cores_ct = 1
                          ) {
  #sparse matrix---------------
  out.all <- lapply(names(ccc_list), function(n) {
    print(n)
    ct_tab <- ccc_list[[n]]
    n.split <- unlist(strsplit(n, "|", fixed = T))
    #cl_list.tmp <- cl_list[[n]]
    out_split <- lapply(n.split, function(x) {
      cl_list_sub <- cl_list[[x]]
      idx <- which(ct_tab[1,] == x)
      cl_list_sub[!names(cl_list_sub) %in% ct_tab[, idx+1]] <- 0
      g.not <- ct_tab[ct_tab$u12 == 0, idx +1]
      g.not <- unique(g.not)
      cl_list_sub[as.character(g.not)] <- 0
      ptw_list_sub <- gs_list[[x]]
      ptw_list_sub[[x]] <- cl_list_sub
      ict_pct_out <- single_integrated_cross_talk(S_list = ptw_list_sub, ref = x,
                             perm_link = perm_link, perm_weights = perm_weights, cut_par = cut_par,
                             A = A, k = k, mc_cores_perm = mc_cores_perm, mc_cores_ct = mc_cores_ct)#complete
      if(length(ict_pct_out) > 1) {
        ict_pct_out <- cbind(ccc = rep(n, nrow(ict_pct_out)), ict_pct_out)
      }
      
      return(ict_pct_out)
    })
    out_split <- out_split[lengths(out_split) > 1]
    out_split <- do.call(rbind, out_split)
    #names(out_split) <- n.split
    return(out_split)
    
    
  })
  out.all <- do.call(rbind, out.all)
  #names(out.all) <- names(ccc_list)
  return(out.all)
  
}
