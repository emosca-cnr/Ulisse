#' Function to obtain DEG ist for multicl_communication()
#' @description `multi_DEGlist()` is used to organize DEG results in a format compatible with `multicl_communication()`
#' @param DEG_res a list of tables with the results of Differentially Expressed Gene (DEG) analysis performed on 
#' each cluster pair. Each element in the list should be named as "cl1|cl2", where cl1 is the one used as numerator 
#' in the fold change calculation in the pair
#' @param log2FCname name of the column with log2 Fold Change value in the DEG result tables
#' @param pval_name name of the column with p-value  in the DEG result tables
#' @param gene_universe vector with the names of the gene of interest. It is used to filter the output and return a
#' DEG list with only the genes under study
#' @param threshold_l2FC threshold used to filter the log2FC. For each cluster pair, genes with a log2FC higher 
#' than this value are considered as highly expressed in the first cluster, otherwise genes lower 
#' than -`threshold_l2FC` are considered as highly expressed in the second cluster
#' @param threshold_pval threshold used in each cluster pair to filter the p-value to consider the genes as significative
#' @details `multi_DEGlist()` takes as inputs data derived from Differentially Expressed Genes (DEG) analysis performed on 
#' each cluster pair and reformat them in a format compatible with `multicl_communication()`. For each cluster pair,
#' the function calls `single_DEGlist` to filter the genes using `threshold_l2FC` and `threshold_pval` and to return
#' a list of two named vector, one for each cluster in the pair. The vectors are the -log10 of the p-value, 
#' named by their gene names and filtered to consider only the one present in `gene_universe`. 
#' The vectors of the clusters are named according to the names of the DEG tables in `DEG_res`.
#' @return The output is composed by a list of objects, one for each cluster pair considered. Each object is named 
#' accordingly to the names in `DEG_res`. Each object is a list of weights named by their corresponding gene names.
#' 

multi_DEGlist <- function(DEG_res=DEG_list, log2FCname = "avg_log2FC", pval_name = "p_val", 
                          threshold_l2FC = 0.25, threshold_pval = 0.1, gene_universe) {
  
  out <- list()
  for (i in 1:length(DEG_res)) {
    n <- unlist(strsplit(names(DEG_res[i]), split = "|", fixed = T))
    tmp <- single_DEGlist(cl = n, gene = rownames(DEG_res[[i]]), 
                          l2FC = DEG_res[[i]][,log2FCname],
                          pval = DEG_res[[i]][,pval_name], gene_universe = gene_universe, threshold_l2FC, threshold_pval)
    out[[i]] <- tmp
    names(out)[i] <- names(DEG_res)[i]
  }
  return(out)
  
}

