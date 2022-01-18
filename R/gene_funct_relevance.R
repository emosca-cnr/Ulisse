#' Function to calculate gene functional relevance
#' @description gene_funct_relevance() calculates functional diversity and interactor diversity for the gene 
#' that participate in significative pathway cross-talks
#' @details gene_funct_relevance() takes as an input the data.frame obtained from `pathway_cross_talk()` and the 
#' adjacency matrix used as an input of `pathway_cross_talk()`. The function then extract the genes partecipating in
#' significative cross-talks and calculate for each of them the functional diversity and the interactor diversity. 
#' For each gene, the functional diversity is the number of pathway with which the gene participate in the formation
#' of a cross-talk of interest; the interactor diversity is the number of different genes with which the gene has links 
#' that contribute to the formation of a cross-talk
#' @param pct the data.frame obtained as a result of `pathway_corss_talk()`
#' @param adj the adjacency matrix used as an input in `pathway_corss_talk()`
#' @param to_plot logical, indicating whether to save the plot of functional relevance
#' @param file_name file name of the plot, if NULL the plot will be saved in the working 
#' directory as "gene_functional_relevance.jpeg"
#' @param plot_names logical, indicating whether to print gene names in the functional relevance plot
#' @return The function returns a data.frame with:
#' \itemize{
#' \item gene: the gene analysed
#' \item nPTW: functional diversity
#' \item nInteractors: interactor diversity
#' \item pathway: list of names of the pathway counted in `nPTW` separated by `;`
#' \item interactors: list of names of the genes counted in `nInteractors` spearated by `;`
#' }
#' @examples 
#' ptw_list <- list(ptwA = c("A", "B","C"), ptwB = c("D", "E", "F"), ptwC = c("A", "B", "E"))
#' adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#' ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#' wgt <- rep(1, 6)
#' pct <- pathway_cross_talk(pathway_list = ptw_list, gene_network_adj = adj, weight = wgt, 
#'   mc_cores_pct = 1, mc_cores_perm = 1, k = 9)
#' funct_rel <- gene_funct_relevance(pct, adj)
#' @import plotrix
#' @export



gene_funct_relevance <- function(pct, adj, to_plot = T, file_name = NULL, plot_names =T) {
  
  ptw.sign <- unique(c(pct$pathway_1, pct$pathway_2))
  genes <- unique(c(unlist(strsplit(pct$gene_pathway1, ";", fixed = T)),
                    unlist(strsplit(pct$gene_pathway2, ";", fixed = T))))
  gene.data <- data.frame(matrix(data=0, nrow = length(genes), ncol = 5, 
                                 dimnames = list(NULL, c("gene", "nPTW", "nInteractors", "pathway", "interactors"))),
                          stringsAsFactors = F)
  gene.data$gene <- genes
  adj <- sign(adj)
  #adjF <- adj[which(rownames(adj) %in% as.character(genes)), which(colnames(adj) %in% as.character(genes))]
  for(i in 1:length(genes)) {
    idx_1 <- grep(genes[i], pct$gene_pathway1)
    idx_2 <- grep(genes[i], pct$gene_pathway2)
    nPTW <- length(unique(c(pct$pathway_2[idx_1], pct$pathway_1[idx_2])))
    ptw <- paste(unique(c(pct$pathway_2[idx_1], pct$pathway_1[idx_2])), collapse = ";")
    
    gene.data[i,2] <- nPTW
    gene.data[i,4] <- ptw
    shr <- unique(c(unlist(strsplit(pct$gene_pathway1[idx_1], ";", fixed = T)),
                    unlist(strsplit(pct$gene_pathway2[idx_2], ";", fixed = T))))
    wlnk <- unique(c(unlist(strsplit(pct$gene_pathway1[idx_2], ";", fixed = T)),
                     unlist(strsplit(pct$gene_pathway2[idx_1], ";", fixed = T))))
    idx <- which(wlnk %in% shr)
    if(length(idx) > 0) {
      wlnk <- wlnk[-idx]
    }
    nInter <- adj[which(rownames(adj) == genes[i]), which(colnames(adj) %in% wlnk), drop = F]
    gene.data[i, 5] <- paste(colnames(nInter)[nInter>0], collapse = ";")
    nInter <- sum(nInter)
    gene.data[i,3] <- nInter
  }
  gene.data <- gene.data[which(gene.data$nPTW > 0),]
  if(to_plot == T) {
    if(is.null(file_name)) {
      if(plot_names == T) {
        grDevices::jpeg("gene_functional_relevance.jpeg", width = 200, height = 200,
             res = 300, units = "mm")
        plot(jitter(gene.data$nPTW, factor = 0.5), jitter(gene.data$nInteractors, factor = 0.5),
             pch=20)
        plotrix::thigmophobe.labels(gene.data$nPTW,
                                    gene.data$nInteractors, 
                                    gene.data$gene, cex=0.5)
        grDevices::dev.off()
      } else {
        grDevices::jpeg("gene_functional_relevance.jpeg", width = 200, height = 200,
             res = 300, units = "mm")
        plot(jitter(gene.data$nPTW, factor = 0.5), jitter(gene.data$nInteractors, factor = 0.5),
             pch=20)
        grDevices::dev.off()
      }
      
    } else {
      if(plot_names == T) {
        grDevices::jpeg(file_name, width = 200, height = 200,
             res = 300, units = "mm")
        plot(jitter(gene.data$nPTW, factor = 0.5), jitter(gene.data$nInteractors, factor = 0.5),
             pch=20)
        plotrix::thigmophobe.labels(gene.data$nPTW,
                                    gene.data$nInteractors, 
                                    gene.data$gene, cex=0.5)
        grDevices::dev.off()
      } else {
        grDevices::jpeg(file_name, width = 200, height = 200,
             res = 300, units = "mm")
        plot(jitter(gene.data$nPTW, factor = 0.5), jitter(gene.data$nInteractors, factor = 0.5),
             pch=20)
        grDevices::dev.off()
      }
    }

  } else {
    NULL
  }

  return(gene.data)
  
}
