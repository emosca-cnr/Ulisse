#' Function to calculate gene functional relevance
#' @description gene_funct_relevance() calculates functional diversity and interactor diversity for the gene 
#' that participate in significative pathway cross-talks
#' @details gene_funct_relevance() takes as an input the data.frame obtained from `pathway_cross_talk()` and the 
#' adjacency matrix used as an input of `pathway_cross_talk()`. The function then extract the genes partecipating in
#' significative cross-talks and calculate for each of them the functional diversity and the interactor diversity. 
#' For each gene, the functional diversity is the number of pathway with which the gene participate in the formation
#' of a cross-talk of interest; the interactor diversity is the number of different genes with which the gene has links 
#' that contribute to the formation of a cross-talk
#' @param pct the data.frame obtained as a result of `pathway_cross_talk()`
#' @param adj the adjacency matrix used as an input in `pathway_cross_talk()`
#' @param only_diff logical, indicating if the calculation of the interactor diversity should include only the genes 
#' with which the gene of interest never share a pathway
#' @param to_plot logical, indicating whether to save the plot of functional relevance
#' @param file_name file name of the plot, if NULL the plot will be saved in the working 
#' directory as "gene_functional_relevance.jpeg"
#' @param plot_names logical, indicating whether to print gene names in the functional relevance plot
#' @param n_cores number of cores to be use to parallelize gene functional relevance analysis
#' @return The function returns a data.frame with:
#' \itemize{
#' \item gene: the gene analysed
#' \item nPTW: functional diversity
#' \item nInteractors: interactor diversity
#' \item pathway: list of names of the pathway counted in `nPTW` separated by `;`
#' \item interactors: list of names of the genes counted in `nInteractors` spearated by `;`
#' }
#' @examples 
#' \dontrun{
#' ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#' adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#' ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#' wgt <- rep(1, 6)
#' pct <- pathway_cross_talk(pathway_list = ptw_list, gene_network_adj = adj, weight = wgt, 
#' genes = LETTERS[1:6], mc_cores_pct = 1, mc_cores_perm = 1, k = 9)
#' funct_rel <- gene_funct_relevance(pct, adj, to_plot=FALSE)
#' }
#' @import plotrix
#' @importFrom grDevices dev.off jpeg
#' @export



gene_funct_relevance <- function(pct, adj, only_diff = T, to_plot = T, 
                                 file_name = NULL, plot_names =T, n_cores = 1) {
  
  ptw.sign <- unique(c(pct$pathway_1, pct$pathway_2))
  gene1 <- unlist(strsplit(pct$gene_pathway1, ";", fixed = T))
  gene2 <- unlist(strsplit(pct$gene_pathway2, ";", fixed = T))
  genes <- unique(c(gene1, gene2))
  genes <- as.character(genes)
  genes <- genes[genes %in% rownames(adj)]
  
  adj <- sign(adj)
  adj <- as.matrix(adj[genes, genes, drop = F])
  
  gene.data <- mclapply(genes, function(i) {
    idx_1 <- grep(i, pct$gene_pathway1)
    idx_2 <- grep(i, pct$gene_pathway2)
    nPTW <- length(unique(c(pct$pathway_2[idx_1], pct$pathway_1[idx_2])))
    ptw <- paste(unique(c(pct$pathway_2[idx_1], pct$pathway_1[idx_2])), collapse = ";")
    wlnk <- unique(c(unlist(strsplit(pct$gene_pathway1[idx_2], ";", fixed = T)),
                     unlist(strsplit(pct$gene_pathway2[idx_1], ";", fixed = T))))
    if(only_diff) {
      shr <- unique(c(unlist(strsplit(pct$gene_pathway1[idx_1], ";", fixed = T)),
                      unlist(strsplit(pct$gene_pathway2[idx_2], ";", fixed = T))))
      
      idx <- which(wlnk %in% shr)
      if(length(idx) > 0) {
        wlnk <- wlnk[-idx]
      }
    }
   
    nInter <- adj[i, colnames(adj) %in% wlnk, drop = F]
    interactors <- paste(colnames(nInter)[nInter>0], collapse = ";")
    nInter <- sum(nInter)
    tmp <- matrix(c(i, nPTW, nInter, ptw, interactors), nrow = 1)
    return(tmp)
  }, mc.cores = n_cores)
  gene.data <- do.call(rbind, gene.data)
  colnames(gene.data) <- c("gene", "nPTW", "nInteractors", "pathway", "interactors")
  gene.data <- data.frame(gene = gene.data[, 1],
                          nPTW = as.numeric(gene.data[, 2]), 
                          nInteractors = as.numeric(gene.data[, 3]),
                          pathway = gene.data[, 4],
                          interactors = gene.data[, 5],
                          stringsAsFactors = F)
  
  gene.data <- gene.data[gene.data$nPTW > 0,]
  if(to_plot == T) {
    if(is.null(file_name)) {
      file_name <- "gene_functional_relevance.jpeg"
      
    } 
    
    if(plot_names == T) {
      jpeg(file_name, width = 200, height = 200,
           res = 300, units = "mm")
      plot(jitter(as.numeric(gene.data$nPTW), factor = 0.5), 
           jitter(as.numeric(gene.data$nInteractors), factor = 0.5),
           pch=20,
           xlab = "Functional diversity", ylab ="Interactor diversity")
      plotrix::thigmophobe.labels(gene.data$nPTW,
                                  gene.data$nInteractors, 
                                  gene.data$gene, cex=0.5)
      dev.off()
    } else {
      jpeg(file_name, width = 200, height = 200,
           res = 300, units = "mm")
      plot(jitter(as.numeric(gene.data$nPTW), factor = 0.5), 
           jitter(as.numeric(gene.data$nInteractors), factor = 0.5),
           pch=20,
           xlab = "Functional diversity", ylab ="Interactor diversity")
      dev.off()
    }
    
  } else {
    NULL
  }
  
  return(gene.data)
  
}
