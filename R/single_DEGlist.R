#' Function to obtain DEG ist for cluster_communication()
#' @description `single_DEGlist()` is used to organize DEG results in a format compatible with `cluster_communication()`
#' @param cl names of the two cluster, should be ordered such that the first cluster is the numerator in the 
#' fold change calculation
#' @param gene vector with the names of the gene of interest
#' @param l2FC log2fold change of the genes in `gene` calculated between the clusters in `cl`
#' @param pval p-value of the genes in `gene` calcualted between the clusters in `cl`
#' @param gene_universe vector with the names of the gene of interest. It is used to filter the output and return a
#' DEG list with only the genes under study
#' @param threshold_l2FC threshold used to filter the log2FC. Genes with a log2FC higher than this value are considered as
#' highly expressed in the first cluster, otherwise genes lower than -`threshold_l2FC` are considered as highly expressed
#' in the second cluster
#' @param threshold_pval threshold used to filter the p-value to consider the genes as significative
#' @details `single_DEGlist()` takes as inputs data derived from Differentially Expressed Genes (DEG) analysis and 
#' reformat them in a format compatible with `cluster_communication()`. In detail, the function considers only the gene
#' with an absolute value of `l2FC` equal or higher than `threshold_l2FC` and a p-value equal or lower than `threshold_pval`.
#' The genes with `l2FC` equal or higher than `threshold_l2FC` are assigned to the first cluster, the ones with 
#' `l2FC` equal or lower than -`threshold_l2FC` are assigned to the second cluster. The function then creates a 
#' named vector for each cluster composed by the -log10 of `pval` of the genes assigned to that clusters, named by their 
#' gene name. The named vectors are then filtered to maintain only the genes present in `gene_universe`.
#'The function returns a list of the two vectors, named with the name in `cl` (in the same order)
#' @return The output is composed by two named vectors, one for each cluster. The vectors are equal to the -log10 of the 
#' p-value of the genes assigned to that cluster, named by the gene names
#' @examples 
#' cl <- c("cl1", "cl2")
#' gene <- letters[1:6]
#' l2FC <- c(0.30, 0.30, 0.30, -0.30, -0.30, -0.30)
#' pval <- rep(0.1, 6)
#' threshold_l2FC <- 0.25
#' threshold_pval <- 0.1
#' universe <- letters[1:6]
#' exp <- single_DEGlist(cl = cl, gene = gene, l2FC = l2FC, pval = pval, 
#'                       threshold_l2FC = 0.25, threshold_pval = 0.1,
#'                       gene_universe = universe)

 


single_DEGlist <- function(cl = c("0", "1"), gene = rownames(DEG_list[[1]]), l2FC = DEG_list[[1]]$avg_log2FC, 
                           gene_universe = unique(c(ligand, receptor)),
                           pval = DEG_list[[1]]$p_val, threshold_l2FC = 0.25, threshold_pval = 0.1) {
  tmp <- data.frame(gene = gene, l2FC = l2FC, pval = pval, stringsAsFactors = F)
  cl1_idx <- which(tmp$l2FC >= threshold_l2FC & tmp$pval <= threshold_pval)
  cl1 <- tmp$pval[cl1_idx]
  names(cl1) <- tmp$gene[cl1_idx]
  cl1 <- -log10(cl1)
  cl1 <- cl1[which(names(cl1) %in% gene_universe)]
  
  cl2_idx <- which(tmp$l2FC <= -(threshold_l2FC) & tmp$pval <= threshold_pval)
  cl2 <- tmp$pval[cl2_idx]
  names(cl2) <- tmp$gene[cl2_idx]
  cl2 <- -log10(cl2)
  cl2 <- cl2[which(names(cl2) %in% gene_universe)]
  
  out <- list(cl1 = cl1, cl2 = cl2)
  names(out) <- cl
  return(out)
}
