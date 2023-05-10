#' Function to calculate gene functional relevance
#' @description The function calculates functional diversity and interactor diversity involved in provided cross-talks
#' @details The function takes as an input the data.frame obtained from `gs_cross_talk()` and the 
#' adjacency matrix used as an input of `gs_cross_talk()`. These inputs are used to obtain the genes involved in
#' significant cross-talks and calculate for each of them the functional diversity and the interactor diversity. 
#' For each gene, the functional diversity is the number of gene-sets (GS) with which the gene is involved in the formation
#' of a cross-talk; the interactor diversity is the number of different genes with which the gene has links 
#' that contribute to the formation of a cross-talk. The function calculate the two measures on the full `gs_cross_talk()` result provided
#' @param ct the data.frame obtained as a result of `gs_cross_talk()`, filtered if needed
#' @param adj the adjacency matrix used as an input in `gs_cross_talk()`
#' @param method = c("count", "relative") If "count" then the number of genes and gene-sets are listed, if "relative" the counts are 
#' relative by using a general model. In this case `ct_null` is required
#' @param ct_null required only if `method = "relative"`. This should be a cross-talk result table coming from a general model
#' @param n_cores number of cores to be use to parallelize gene functional relevance analysis
#' @return If `method = "count` The function returns a data.frame with:
#' \itemize{
#' \item gene: the gene analysed
#' \item functional_diversity: functional diversity
#' \item interactor_diversity: interactor diversity
#' \item functional relevance: log2(`functional_diversity`/`interactor_diversity`)
#' \item n_gs_gene: number of gs in which `gene` is present
#' \item gs_gene: GS in which `gene` is present
#' \item functional_gs: list of names of the gene-sets counted in `functional_diversity` separated by `;`
#' \item interactors_gene: list of names of the genes counted in `interactor_diversity` separated by `;`
#' } 
#' Otherwise, if `method = "relative` the function returns a data.frame with:
#' \itemize{
#' \item gene: the gene analysed
#' \item relative_functional_diversity: relative functional diversity
#' \item relative_interactor_diversity: relative interactor diversity
#' \item relative_functional relevance: log2(`relative_functional_diversity`/`relative_interactor_diversity`)
#' \item functional_diversity: functional diversity calculated from `ct`
#' \item interactor_diversity: interactor diversity calculated from `ct`
#' \item functional relevance: log2(`functional_diversity`/`interactor_diversity`)
#' \item n_gs_gene: number of gs in which `gene` is present
#' \item functional_gs: list of names of the gene-sets counted in `functional_diversity` separated by `;`
#' \item interactors_gene: list of names of the genes counted in `interactor_diversity` separated by `;`
#' \item functional_diversity_null: functional diversity calculated from general model `ct_null`
#' \item interactor_diversity_null: interactor diversity calculated from general model `ct_null`
#' \item functional_gs_null: list of names of the gs counted in `functional_diversity_null` separated by `;`
#' \item interactors_gene_null: list of names of the genes counted in `interactor_diversity_null` separated by `;`
#' }
#' @examples 
#' \dontrun{
#' gslist <- list(gsA = c("A", "B","C"), gsB = c("D", "E", "F"), gsC = c("A", "B", "E"))
#' adj <- matrix(data = sample(c(0,1), 6*6, replace = TRUE), nrow = 6, 
#' ncol = 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#' wgt <- rep(1, 6)
#' pct <- gscross_talk(gslist = gslist, gene_network_adj = adj, mc_cores_pct = 1, mc_cores_perm = 1, k = 9)
#' funct_rel <- gene_funct_relevance(ct, adj, to_plot=FALSE)
#' }
#' @import plotrix
#' @importFrom grDevices dev.off jpeg
#' @export



gene_functional_relevance <- function(ct, adj, method = "count", ct_null = NULL, n_cores = 1) {
  
  ptw.sign <- unique(c(ct$gs1, ct$gs2))
  gene1 <- unlist(strsplit(ct$genes_gs1, ";", fixed = T))
  gene2 <- unlist(strsplit(ct$genes_gs2, ";", fixed = T))
  genes <- unique(c(gene1, gene2))
  genes <- as.character(genes)
  genes <- genes[genes %in% rownames(adj)]
  adjOrg <- adj
  adj <- sign(adj)
  adj <- as.matrix(adj[genes, genes, drop = F])
  
  gene.data <- mclapply(genes, function(i) {
    idx_1 <- grep(i, ct$genes_gs1)
    idx_2 <- grep(i, ct$genes_gs2)
    nCl_in <- length(unique(c(ct$gs1[idx_1], ct$gs2[idx_2])))
    ptw_in <- paste(unique(c(ct$gs1[idx_1], ct$gs2[idx_2])), collapse = ";")
    nPTW <- length(unique(c(ct$gs2[idx_1], ct$gs1[idx_2])))
    ptw <- paste(unique(c(ct$gs2[idx_1], ct$gs1[idx_2])), collapse = ";")
    wlnk <- unique(c(unlist(strsplit(ct$genes_gs1[idx_2], ";", fixed = T)),
                     unlist(strsplit(ct$genes_gs2[idx_1], ";", fixed = T))))
    
   
    nInter <- adj[i, colnames(adj) %in% wlnk, drop = F]
    interactors <- paste(colnames(nInter)[nInter>0], collapse = ";")
    nInter <- sum(nInter)
    tmp <- matrix(c(i,nCl_in, nPTW, nInter, ptw_in, ptw, interactors), nrow = 1)
    return(tmp)
  }, mc.cores = n_cores)
  gene.data <- do.call(rbind, gene.data)
  # colnames(gene.data) <- c("gene", "nGS", "nInteractors", "gs", "interactors")
  gene.data <- data.frame(gene = gene.data[, 1],
                          functional_diversity = as.numeric(gene.data[, 3]), 
                          interactor_diversity = as.numeric(gene.data[, 4]),
                          functional_relevance = log2(as.numeric(gene.data[, 3])/as.numeric(gene.data[, 4])),
                          n_gs_gene = as.numeric(gene.data[, 2]),
                          gs_gene = gene.data[,5],
                          functional_gs = gene.data[, 6],
                          interactors_gene = gene.data[, 7],
                          stringsAsFactors = F)
  #gene.data <- gene.data[gene.data$functional_diversity > 0,]
  
  if(method == "relative") {
    #null model to normalize--------------------------------
    ptw.sign.null <- unique(c(ct_null$gs1, ct_null$gs2))
    gene1.null <- unlist(strsplit(ct_null$genes_gs1, ";", fixed = T))
    gene2.null <- unlist(strsplit(ct_null$genes_gs2, ";", fixed = T))
    genes.null <- unique(c(gene1.null, gene2.null))
    genes.null <- as.character(genes.null)
    genes.null <- genes.null[genes.null %in% rownames(adjOrg)]
    
    adj.null <- sign(adjOrg)
    adj.null <- as.matrix(adj.null[genes.null, genes.null, drop = F])
    
    gene.data.null <- mclapply(genes.null, function(i) {
      idx_1 <- grep(i, ct_null$genes_gs1)
      idx_2 <- grep(i, ct_null$genes_gs2)
      nPTW <- length(unique(c(ct_null$gs2[idx_1], ct_null$gs1[idx_2])))
      ptw <- paste(unique(c(ct_null$gs2[idx_1], ct_null$gs1[idx_2])), collapse = ";")
      wlnk <- unique(c(unlist(strsplit(ct_null$genes_gs1[idx_2], ";", fixed = T)),
                       unlist(strsplit(ct_null$genes_gs2[idx_1], ";", fixed = T))))
        
      nInter <- adj.null[i, colnames(adj.null) %in% wlnk, drop = F]
      interactors <- paste(colnames(nInter)[nInter>0], collapse = ";")
      nInter <- sum(nInter)
      tmp <- matrix(c(i, nPTW, nInter, ptw, interactors), nrow = 1)
      return(tmp)
    }, mc.cores = n_cores)
    gene.data.null <- do.call(rbind, gene.data.null)
    #colnames(gene.data.null) <- c("gene", "nPTW_null", "nInteractors_null", "pathway_null", "interactors_null")
    gene.data.null <- data.frame(gene = gene.data.null[, 1],
                                 functional_diversity_null = as.numeric(gene.data.null[, 2]), 
                                 interactor_diversity_null = as.numeric(gene.data.null[, 3]),
                                 functional_gs_null = gene.data.null[, 4],
                                 interactors_gene_null = gene.data.null[, 5],
                                 stringsAsFactors = F)
    gene.data <- merge(gene.data, gene.data.null, by = "gene", all.x = T)
    gene.data$relative_functional_diversity <- gene.data$functional_diversity/gene.data$functional_diversity_null
    gene.data$relative_interactor_diversity <- gene.data$interactor_diversity/gene.data$interactor_diversity_null
    gene.data$relative_functional_relevance <- log2(gene.data$relative_functional_diversity/gene.data$relative_interactor_diversity)
    gene.data <- gene.data[, c("gene", "relative_functional_diversity", "relative_interactor_diversity", "relative_functional_relevance",
                               "functional_diversity", "interactor_diversity", "functional_relevance", "n_gs_gene", 
                               "gs_gene", "functional_gs", "interactors_gene",
                               "functional_diversity_null", "interactor_diversity_null",
                               "functional_gs_null", "interactors_gene_null")]
  }
  
  
  
  return(gene.data)
  
}
