#' Function to calculate cross-talks between pathways of different gene communities
#' @description Calculates cross-talks between pathways composed by gene of different gene communities
#' @details The function uses the `membership` data to re-organize the `pathway_list` as to obtain a list of 
#' pathway list for each gene community. To do so genes will be grouped into communities and then grouped into pathways.
#' At this point, the function calculates CT between each pathway pairs of different 
#' communities that shows at least a link.
#' @param gs_list a named list of genes grouped into pathways
#' @param gene_network_adj gene network adjacency matrix
#' @param weight weights of the genes in `genes`. If not provided, the function assigns to each gene
#' a weight of 1
#' @param membership named vector with the membership of `genes` to topological communities. If set to `NULL`
#' the function calculates the communities by using fastgreedy algorithm
#' @param genes vector with the genes of interest to be used for TM-PCT calculation
#' @param mc_cores_ct number of threads to be used to calculate cross talk 
#' @param mc_cores_tm number of threads to be used to calcualte TM-PCT on different communities combination
#' @return The function returns a list of two object:
#' \enumerate{
#'  \item comm_pathway_list: the communities pathway lists used for TM-PCT calculation. Each list is named with the name 
#'  of the gene community
#'  \item {TM_PCT_res}: a data frame with
#' \itemize{
#'  \item commID_1, commID_2: community to which belongs the genes in `pathway_1` and `pathway_2`, respectively
#'  \item pathway_1, pahtway_2: the pathway pair considered
#'  \item score: TM-PCT score
#'  \item ngenes_pathway_1, ngenes_pathway_2: number of genes involved in `pathway_1` and `pathway_2`, respectively
#'  \item n_link: number of links between the pathways considered
#'  \item weight_pathway1, weight_pathway2: cumulative weights of the genes involved in CT
#'  in `pathway_1` and `pathway_2`
#'  \item gene_pathway1, gene_pathway2: gene involved in the CT in `pathway_1` and `pathway_2`, respectively
#'  }
#'  }
#'  If membership = NULL, then also the membership is returned
#' @examples  
#' \dontrun{
#'  ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#'  adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#'  ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#'  wgt <- rep(1, 6)
#'  memb <- c(1, 1, 2, 2, 3, 3)
#'  pct <- TM_PCT(pathway_list = ptw_list, gene_network_adj = adj, weight = wgt, membership = memb, 
#'                 mc_cores_tm = 1, mc_cores_ct = 1)
#' }
#' @import parallel
#' @importFrom gtools permutations
#' @import kit
#' @importFrom stringi stri_c
#' @importFrom methods is as
#' @export

TM_CT <- function (gs_list, gene_network_adj, membership, genes, weight, 
                    mc_cores_ct = 2, mc_cores_tm = 2) {
  
  if(!is(gene_network_adj, "sparseMatrix" )) {
    gene_network_adj <- as(gene_network_adj, "dgCMatrix")
  }
  gene_network_adj <- sign(gene_network_adj)
  gs_list <- lapply(gs_list, function(x) x <- x[names(x) %in% rownames(gene_network_adj)])
  gene <- as.character(unique(unlist(lapply(gs_list, names))))
  #sub_adj_mt <- gene_network_adj[gene, gene, drop = F]
  
  if(is.null(membership)) {
    sub.g <- graph_from_adjacency_matrix(gene_network_adj, mode = "undirected")
    community <- fastgreedy.community(sub.g)
    membership <- membership(community)
    ret_memb <- T
    
  } else {
    ret_memb <- F
  }
   
  unqMem <- unique(membership)
  mem_ptwL <- mclapply(unqMem, function(z) {
    tmp <- lapply(gs_list, function(x) x <- x[names(x) %in% rownames(gene_network_adj)[membership == z]])
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
  }, mc.cores = mc_cores_ct
  )
  idx <- which(mclapply(xx, function(n) sum(n, na.rm = T), mc.cores = mc_cores_ct) !=0)
  if(length(idx) == 0) {
    print("No cambination of TM available to calculate PCT") 
    return("No cambination of TM available to calculate PCT")
  } else {
    xx <- xx[idx]
    comb_mem <- comb_mem[idx,,drop = F]
    
    
    res <- parallel::mclapply(1:nrow(comb_mem), function(j) {
      ptw_listCC1 <- mem_ptwL[[comb_mem[j,1]]]
      ptw_listCC2 <- mem_ptwL[[comb_mem[j,2]]]
      
      comb_p <- expand.grid(ptwCC1 = names(ptw_listCC1),
                            ptwCC2 = names(ptw_listCC2), stringsAsFactors = F)
      comb_p <- comb_p[comb_p$ptwCC1 != comb_p$ptwCC2,]
      
      xxCC <- mclapply(1:nrow(comb_p),  function (x) {
        idx_1 <- as.character(names(ptw_listCC1[[comb_p[x,1]]]))
        idx_2 <- as.character(names(ptw_listCC2[[comb_p[x,2]]]))
        tmp <- gene_network_adj[idx_1, idx_2, drop = F]
        
        return( tmp) 
      }, mc.cores = mc_cores_ct)
      
      
      idxCC <- which(mclapply(xxCC, function(n) sum(n, na.rm = T), mc.cores = mc_cores_ct) !=0)
      if (length(idxCC) != 0) {
        xxCC <- xxCC[idxCC]
        comb_p <- comb_p[idxCC,]
        
        pct <- mclapply(1:length(xxCC), function(m) {
          tab <- xxCC[[m]]
          
          g.1 <- ptw_listCC1[[comb_p[m, 1]]]
          g.2 <- ptw_listCC2[[comb_p[m, 2]]]
          
          tab.out <- cross_talk(mat = tab, weight = list(g1 = g.1, g2 = g.2))
          
          tab.out <- matrix(c(comb_mem[j, 1], comb_p[m,1], comb_mem[j, 2], comb_p[m,2], tab.out), nrow = 1)
          
          return(tab.out)
        },   mc.cores = mc_cores_ct)
        pct <- do.call(rbind, pct)
        # colnames(pct) <- c("commID_1", "gs1", "commID_2", "gs2", "ct_score", 
        #                    "ngenes_pathway1", "ngenes_pathway2", "nlink", "weight_pathway1", 
        #                    "weight_pathway2", "gene_pathway1", "gene_pathway2")
        
        return(pct)
      } else {
        ans <- NULL
        return(ans)
      }
      
    }, mc.cores = mc_cores_tm)
    
    res <- do.call(rbind, res)
    res <- data.frame(commID_1 = res[, 1],
                      gs1 = res[, 2], 
                      commID_2 = res[, 3], 
                      gs2 = res[, 4], 
                      ct_score = as.numeric(res[, 5]), 
                      ngenes_pathway1 = as.numeric(res[, 6]), 
                      ngenes_pathway2 = as.numeric(res[, 7]), 
                      nlink = as.numeric(res[, 8]), 
                      weight_pathway1 = as.numeric(res[, 9]), 
                      weight_pathway2 = as.numeric(res[, 10]), 
                      gene_pathway1 = res[, 11], 
                      gene_pathway2 = res[, 12], stringsAsFactors = F)
    
    if(ret_memb) {
      out <- list(membership = membership ,comm_pathway_list = mem_ptwL, TM_CT_res = res)
      return(out)
    } else {
      out <- list(comm_pathway_list = mem_ptwL, TM_CT_res = res)
      return(out)
    }
    
  }
  
}
