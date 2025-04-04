#' Integrated cross talk on specific cluster analysis
#' @description The function calculates cross-talk between genes involved in communication with a specific cluster/cell type and cluster gene sets
#' @details The function takes as inputs the adjacency matrix of the biological network (`A`), the gene-set list `S_list()` and the name of the gene set that
#'  contains the communication genes (`ref`). For each gene-set pair that shows at least a link between them, the function calculates the cross-talk 
#'  and the statistical evaluation, as described in the paper, using `k` permutations. 
#' @param S_list a named list of genes grouped into gene sets, as obtained from `build_S_list()`, plus the reference gene set list
#' @param ref name of the gene set in `S_list` that contains all genes of cluster i involved in the communication with cluster j
#' @param A adjacency matrix of the whole gene network considered (can be a sparseMatrix)
#' @param k number of permutations
#' @param bin_type can be either "number" (suggested), "interval" or "width". See `ggplot2::cut_interval()` for details.
#' @param cut_par number of bin to cut. If set to `NULL`, the function will search for the best cut in 2:15 
#' @param perm_link = c("degree", "simple") If the permutation of the link should be degree-conservative ("degree"), or random ("simple")
#' @param perm_weights = c("degree", "simple") If the permutation of the weights should be degree-conservative ("degree"), or random ("simple")
#' @param mc_cores_ct number of threads to be used for cross talk calculation
#' @param mc_cores_perm number of thread to be used in permutations
#' @return The function returns a data.frame with the results of ct calculation:
#' \itemize{
#'    \item S1_name,S2_name: names of gene-set pairs considered
#'    \item c: cross-talk score
#'    \item S1_size,S2_size: number of genes in S1 and S2
#'    \item S1_s2_size,S2_s1_size: number of S1 genes interacting with S2, and vice versa
#'    \item dL: number of links between S1 and S2
#'    \item L: number of possible links between S1 and S2
#'    \item r_c: cross-talk saturation, calculated as `dL/L`
#'    \item u1,u2: sum of the gene weights in S1 and S2, respectively
#'    \item S1,S2: list of interacting genes in S1 and S2, respectively
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
#' @importFrom NPATools perm_vertices perm_X0 calc_p
#' @export

single_integrated_cross_talk <- function(S_list, ref, A, k = 9,
                    bin_type = "number", cut_par = 9, perm_link = "degree", #c("degree", "simple")
                    perm_weights = "simple", #c("degree", "simple")
                    mc_cores_perm = 1, mc_cores_ct = 1) {
  #initial check---------------
  perm_link <- match.arg(perm_link, c("simple", "degree"))
  perm_weights <- match.arg(perm_weights, c("simple", "degree"))
  bin_type <- match.arg(bin_type, c("interval", "number", "width"))
  #ct_type <- match.arg(ct_type, c("pct", "ccc"))
  S_list <- lapply(S_list, function(x) x <- x[names(x) %in% rownames(A)])
  #dense to sparse M------------------
  A <- sign(as.matrix(A))
  
  #finding cross-talk > 0-----------------
  n.gs <- names(S_list)
  n.gs <- n.gs[!n.gs %in% ref]
  comb_p <- expand.grid(ref, names(S_list), stringsAsFactors = F)
  
  cat("#{Gene set pairs}:", nrow(comb_p), "\n")
  
  #Adjancency for each pair S_i, S_j
  sub_adj <- mclapply(1:nrow(comb_p), function(x) {
    g.1 <- get_elem(S_list, comb_p[x, 1])
    n.g1 <- as.character(names(g.1))
    g.2 <- get_elem(S_list, comb_p[x, 2])
    n.g2 <- as.character(names(g.2))
    
    s <- intersect(n.g1, n.g2)
    
    g.1 <- g.1[!g.1 == 0]
    g.2 <- g.2[!g.2 == 0]
    n.g1 <- as.character(names(g.1))
    n.g2 <- as.character(names(g.2))
    
    n.g1 <- n.g1[!n.g1 %in% s]
    n.g2 <- n.g2[!n.g2 %in% s]
    
    tab <- ss(A, n.g1, n.g2)
    tab <- sum(tab, na.rm = T)  
    return(tab)
  }, mc.cores = mc_cores_ct)
  
  #removes null A
  idx <- which(unlist(sub_adj) > 0)
  if(length(idx) == 0) {
    print("no available CT")
    return("no available CT")
  } else {
    comb_p <- comb_p[idx,, drop = F]
    ptw <- unique(c(comb_p[, 1], comb_p[, 2]))
    S_list <- S_list[ptw]
    
    #using only genes in S_list------------------
    g <-  unique(as.character(names(unlist(unname(S_list)))))
    A <- ss(A, g, g)
    rm(g)
    
    #### REAL CT
    cat("Calculation of cross-talks...\n")
    ans <- mclapply(1:nrow(comb_p), function(z){
      
      #gene set names
      n.gs1 <- comb_p[z, 1]
      n.gs2 <- comb_p[z, 2]
      
      #gene set genes
      gs.1w <- get_elem(S_list, n.gs1)
      gs.2w <- get_elem(S_list, n.gs2)
      
      n.g1 <- names(gs.1w)
      n.g2 <- names(gs.2w)
      
      #remove the intersection
      s <- intersect(n.g1, n.g2)
      n.g1 <- n.g1[!n.g1 %in% s]
      n.g2 <- n.g2[!n.g2 %in% s]
      gs.1w <- gs.1w[n.g1]
      gs.2w <- gs.2w[n.g2]
      
      #links
      tab.w <- ss(A, n.g1, n.g2)
      
      ct.w <- cross_talk.opt(mat = tab.w, w1 = as.matrix(gs.1w), w2 = as.matrix(gs.2w))
      ct.w <- cbind(S1_name = n.gs1, S2_name = n.gs2, ct.w)
      
      return(ct.w)
      
    })
    
    #real cross-tlak
    ct <- rowbind(ans)
    
    
    #matrix adjusting + perms---------------------
    if(k> 0) {
      cat("Preparing permutations...\n")
      deg.adj <- rowSums(A)
      if(perm_link == "degree" & perm_weights == "degree" | 
         perm_link == "simple" & perm_weights == "degree" | 
         perm_link == "degree" & perm_weights == "simple") {
        if(is.null(cut_par)) {
          cat("cut_par == NULL\nSearching for best cut_par in 2:15...\n")
          n <- c()
          for(j in 2:15) {
            brk <- breaks(x = deg.adj, nbins = j, equal = "number")
            n[j] <- anyDuplicated(brk)
          }
          
          idx <- which(n == 0)
          cut_par <- idx[length(idx)]
          if(length(cut_par) == 0 ) {
            cat("impossible to cut degree in bins. Switching to 'simple' method\n")
            perm_link <- "simple" 
            perm_weights <- "simple"
          } else {
            cat(paste0("cut_par used = ", cut_par, "\n"))
          }
          
          #print(cut_par)
        }
        
      }
      if(perm_link != perm_weights) {
        perm_v <- perm_vertices(vert_deg = deg.adj, k = k, method = perm_weights, 
                                cut_par = cut_par, bin_type = bin_type)
        
        m <- matrix(data = 0, nrow = nrow(A), ncol = length(S_list), 
                    dimnames = list(rownames(A), names(S_list)))
        for(i in 1:length(S_list)) {
          l <- S_list[[i]]
          m[names(l), i] <- l
        }
        perm_m <- perm_X0(X0 = m, perms = perm_v)
        
        perm_v <- perm_vertices(vert_deg = deg.adj, k = k, method = perm_link, 
                                cut_par = cut_par, bin_type = bin_type)
      } else { # permutation over vertices equal to those over links
        perm_v <- perm_vertices(vert_deg = deg.adj, k = k, method = perm_link, 
                                cut_par = cut_par, bin_type = bin_type)
        
        m <- matrix(data = 0, nrow = nrow(A), ncol = length(S_list), 
                    dimnames = list(rownames(A), names(S_list)))
        for(i in 1:length(S_list)) {
          l <- S_list[[i]]
          m[names(l), i] <- l
        }
        perm_m <- perm_X0(X0 = m, perms = perm_v)
        
      }
      rm(deg.adj)
      
      #index of every permutation
      perm_v <- lapply(perm_v, function(l) l <- setNames(1:length(l), l))
      
      
      perm_ans <- lapply(2:length(perm_v) , function(x) { 
        ans <- parallel::mclapply(1:nrow(comb_p), function(z) {
          #print(z)
          #gene set names
          n.gs1 <- comb_p[z, 1]
          n.gs2 <- comb_p[z, 2]
          
          #gene set genes
          gs.1w <- get_elem(S_list, n.gs1)
          gs.2w <- get_elem(S_list, n.gs2)
          
          n.g1 <- names(gs.1w)
          n.g2 <- names(gs.2w)
          
          #remove the intersection
          s <- intersect(n.g1, n.g2)
          n.g1 <- n.g1[!n.g1 %in% s]
          n.g2 <- n.g2[!n.g2 %in% s]
          
          
          gs.1w <- ss(perm_m[[x]], n.g1, n.gs1)
          gs.2w <- ss(perm_m[[x]], n.g2, n.gs2)
          gs.1l <- ss(perm_m[[1]], n.g1, n.gs1)
          gs.2l <- ss(perm_m[[1]], n.g2, n.gs2)
          
          ### permutations over links
          g1.idx.l <- ss(perm_v[[x]], n.g1)
          g2.idx.l <- ss(perm_v[[x]], n.g2)
          
          tab.l <- ss(A, g1.idx.l, g2.idx.l)
          rownames(tab.l) <- n.g1
          colnames(tab.l) <- n.g2
          
          ct.l <- cross_talk.perm(mat = tab.l, w1 = gs.1l, w2 = gs.2l)
          
          ### real links
          tab.w <- ss(A, n.g1, n.g2)
          
          ct.w <- cross_talk.perm(mat = tab.w, w1 = gs.1w, w2 = gs.2w)
          
          return(data.frame(ct.l, ct.w))
          
          
        }, mc.cores = mc_cores_perm)
        
        ans <- rowbind(ans)
        
        return(ans)
      })
      
      perm_ans <- c(list(cbind(ct$c, ct$c)), perm_ans)  
      cat("P-values...\n")
      
      p_val_list_l <- lapply(perm_ans, function(l) l[, 1, drop=F])
      p_val_l <- calc_p(X = p_val_list_l)
      
      p_val_list_w <- lapply(perm_ans, function(l) l[, 2, drop=F])
      
      p_val_w <- calc_p(X = p_val_list_w)
      
      cat("Assembling output...\n")
      ans <- ct
      
      p_comb <- p_val_l * p_val_w
      p_comb <- p_comb - (p_comb *log(p_comb))
      ct_score <- ans$c * -log10(p_comb)
      
    } else{
      
      ct_score <- p_val_l <- p_val_w <- p_comb <- NA
      ans <- ct
    }  
    ans <- cbind(ans, s=as.numeric(ct_score), pA=as.numeric(p_val_l), pU=as.numeric(p_val_w), p=as.numeric(p_comb))
    
    return(ans)
  }
  
}
