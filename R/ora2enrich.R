#' ora2enrich
#' this function translates the result of an ORA into the DOSE class "enrichResult"
#' @param ulisse_res ora analysis result
#' @param pvalueCutoff cut off over adjusted p value= 0.25
#' @param pAdjustMethod method for p value adjustement
#' @param qvalueCutoff cut off over q values
#' @param gene input gene list
#' @param universe universe of genes
#' @param geneSets list of gene sets
#' @param organism organism
#' @param keytype keytype
#' @param ontology ontology
#' @param readable readable
#'
#' @export
#' @importClassesFrom DOSE enrichResult
#' @importFrom qvalue qvalue


ora2enrich <- function(ulisse_res, pvalueCutoff = 0.25, pAdjustMethod = "BH", qvalueCutoff = 0.25, gene = NULL, universe = NULL, geneSets = NULL, organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN", readable = FALSE){

  geneID <- unlist(sapply(ulisse_res$wbd_symbols, function(x) gsub(";", "/", x)))

  if(is.null(ulisse_res$q_val)){
    q_val <- qvalue::qvalue(p=ulisse_res$p, lambda=0.05, pi0.method="bootstrap")$qvalues
  }else{
    q_val <- ulisse_res$q_val
  }

  Over <- data.frame(ID = as.character(ulisse_res$gsid), GeneRatio = paste0(ulisse_res$wbd, "/",  ulisse_res$bd), BgRatio = paste0(ulisse_res$wb, "/",  ulisse_res$N), pvalue = ulisse_res$p, p.adjust=ulisse_res$p_adj, qvalue=q_val, Count=ulisse_res$wbd, geneID=geneID, stringsAsFactors = FALSE)

  Over$Description <- ulisse_res$name
  row.names(Over) <- as.character(Over$ID)

  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff,
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
           gene = as.character(gene), universe = universe, geneSets = geneSets,
           organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN",
           readable = FALSE)

  return(x)


}
