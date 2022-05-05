#' Function to prepare the input for `cluster_communication()` by using DEG analysis results data
#' @description The function prepares the input for `cluster_communication()` analysis by using DEG results data
#' @details The function prepares the gene list in a format compatible with the `cluster_communication()` analysis. 
#' In detail, the function needs some information obtained from DEG analysis: cluster tested, gene names, p-values,
#' log2(Fold Change). These data will be filtered to consider only the gene of interest in the study (`universe`). 
#' The function calculates a score for each gene as: abs(log2FC) * -log10(p-value). 
#' If the p-value is 0, this is approximated to 1/10 of the minimum non-zero p-value present in the data. 
#' The function returns a gene set list for each cluster, each composed by a vector of scores named after the genes 
#' names provided.
#' @param cluster column of the result table of DEG analysis with the data of the cluster under analysis
#' @param p_val column of the result table of DEG analysis with the data of the p-value of the test. Could be either
#' nominal or adjusted p-value
#' @param log2FC column of the result table of DEG analysis with the data of the log2Fold Change
#' @param gene column of the result table of DEG analysis with the names of the genes tested
#' @param universe vector with the names of the gene of interest in teh study. This value should be the names
#' in the LR network
#' @return The function returns a gene list composed by a vector for each cluster provided. Each vector is 
#' composed by the score, calulcated as abs(log2FC) * -log10(p-value), named by the respective genes
#' @export


preparing_DEG_list <- function(cluster, p_val, log2FC, gene, universe) {
  
  out <- data.frame(cluster = cluster, 
                    p_val = p_val, 
                    log2FC = log2FC, 
                    gene = gene, stringsAsFactors = F)
  out <- out[out$gene %in% universe,]
  
  out$p_val[out$p_val == 0] <- min(out$p_val[out$p_val != 0])/10
  
  out$score <- abs(out$log2FC) * -log10(out$p_val)
  list_out <- split(out$score, out$cluster)
  names_list <- split(out$gene, out$cluster)
  list_out <- lapply(1:length(list_out), function(x) {
    tmp <- list_out[[x]]
    names(tmp) <- names_list[[x]]
    return(tmp)
  })
  names(list_out) <- unique(cluster)
  return(list_out)
}

