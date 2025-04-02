#' Cross-talk gene classification analysis
#' @description The function calculates (saturation of) cross-talk diversity and (saturation of) 
#'  interaction diversity involved in provided cross-talks
#' @details The function takes as an input the data.frame obtained from `gs_cross_talk()` and the 
#'  adjacency matrix used as an input of `gs_cross_talk()`. These inputs are used to calculate 
#'  (saturation of) cross-talk diversity and (saturation of) interaction diversity
#'  involved in significant cross-talks, as defined in the paper. 
#'  This function has two applications:
#'  \enumerate{
#'  \item calculate cross-talk diversity and interaction diversity in a reference model
#'  \item calculate saturation of cross-talk diversity and saturation of interaction diversity 
#'    by using the reference obtained in `1`
#'  }
#' @param gs_ct the data.frame obtained as a result of `gs_cross_talk()`, filtered if needed
#' @param A the adjacency matrix used as an input in `gs_cross_talk()`
#' @param ct_ref =`NULL` or a data.frame with cross-talk calculated on a reference. 
#'  If provided, saturation `r_g` and `r_d` are calculated. 
#' @param n_cores number of cores to be use to parallelize gene gene set relevance analysis
#' @return If `ct_ref = NULL` The function returns a data.frame with:
#' \itemize{
#' \item gene: the gene analysed
#' \item d_q: number of altered cross-talk `gene` contributes to
#' \item d_d: number of interactors contributing to any altered cross-talk together with `gene`
#' \item d_q_S: list of names of the gene-sets counted in `d_q` separated by `;`
#' \item d_d_gene: list of names of the genes counted in `d_d` separated by `;`
#' } 
#' Otherwise, if `ct_ref` is provided the function returns a data.frame with:
#' \itemize{
#' \item gene: the gene analysed
#' \item d_q: cross-talk diversity, or number of gene set in altered cross-talk that contains interactors of `gene`
#' \item d_d: interaction diversity, or number of interactors of `gene` that are part of altered cross-talks
#' \item r_q: saturation of cross-talk diversity
#' \item r_d: saturation of interaction diversity
#' \item d_q_S: list of gene set names counted in `d_q` separated by `;`
#' \item d_d_gene: list of names of the genes counted in `d_d` separated by `;`
#' \item q: number of gene set in cross-talk that contains interactors of `gene`
#' \item d: number of interactors of `gene` that are part of cross-talks
#' \item q_S: list of gene set names counted in `q` separated by `;`
#' \item d_gene: list of names of the genes counted in `d` separated by `;`

#' }
#' @importFrom collapse ss rowbind funique
#' @importFrom stringi stri_split
#' @importFrom parallel mclapply
#' @export


gene_classification <- function(gs_ct=NULL, A=NULL, ct_ref = NULL, n_cores = 1) {
  
  
  cat("Preparing data...\n")
  
  #colnames(ct)[match(c("S1_name", "S2_name", "S1", "S2"), colnames(ct))] <- c("gs1", "gs2", "genes_gs1", "genes_gs2") #backward compatibility
  
  ptw.sign <- unique(c(gs_ct$S1_name, gs_ct$S1_name))
  gene1 <- unlist(stri_split(gs_ct$S1, regex = ";"))
  gene2 <- unlist(stri_split(gs_ct$S2, regex = ";"))
  genes <- unique(c(gene1, gene2))
  genes <- as.character(genes)
  genes <- genes[genes %in% rownames(A)]
  
  A <- sign(as.matrix(A))
  A <- ss(A, genes, genes)
  
  cat("Loop through", length(genes), "genes...\n")
  
  gene.data <- mclapply(genes, function(i) {
    
    #cat(".")
    
    idx_1 <- grep(paste0(";", i, ";"), paste0(";", gs_ct$S1, ";")) #; before and after every gene and at the end of ct$genes_gs1
    idx_2 <- grep(paste0(";", i, ";"), paste0(";", gs_ct$S2, ";")) #; before and after every gene  and at the end of ct$genes_gs2
    
    nCl_in <- length(funique(c(gs_ct$S1_name[idx_1], gs_ct$S2_name[idx_2])))
    
    ptw_in <- paste(funique(c(gs_ct$S1_name[idx_1], gs_ct$S2_name[idx_2])), collapse = ";")
    
    nPTW <- length(funique(c(gs_ct$S2_name[idx_1], gs_ct$S1_name[idx_2])))
    
    ptw <- paste(funique(c(gs_ct$S2_name[idx_1], gs_ct$S1_name[idx_2])), collapse = ";")
    
    wlnk <- funique(c(unlist(stri_split(gs_ct$S1[idx_2], regex = ";")),
                     unlist(stri_split(gs_ct$S2[idx_1], regex = ";"))))
    
    nInter <- ss(A, i, colnames(A) %in% wlnk)
    interactors <- paste(colnames(nInter)[nInter>0], collapse = ";")
    nInter <- sum(nInter)
    tmp <- data.frame(i, nPTW, nInter, ptw, interactors, stringsAsFactors = F)
    
    return(tmp)
    
  }, mc.cores = n_cores)
  
  gene.data <- rowbind(gene.data)
  
  colnames(gene.data) <- c("gene", "d_q", "d_d", "d_q_S", "d_d_gene")
  # gene.data <- data.frame(
  #   gene = gene.data[, 1],
  #   d_q = as.numeric(gene.data[, 2]), 
  #   d_d = as.numeric(gene.data[, 3]),
  #   d_q_S = gene.data[, 4],
  #   d_d_gene = gene.data[, 5],
  #   stringsAsFactors = F
  # )
  
  if(!is.null(ct_ref)) {
    
    #ct_ref <- ct_ref[, c(1:3, 6:7)]
    colnames(ct_ref) <- c("gene", "q", "d", "q_S", "d_gene")
    gene.data <- merge(gene.data, ct_ref, by = "gene", all.x = T)
    gene.data$r_q <- gene.data$d_q/gene.data$q
    gene.data$r_d <- gene.data$d_d/gene.data$d
    gene.data <- gene.data[, c("gene","r_q", "r_d",  "d_q", "d_d", "d_q_S", "d_d_gene",
                               "q", "d",
                               "q_S", "d_gene")]
  }
  
  
  
  return(gene.data)
  
}
