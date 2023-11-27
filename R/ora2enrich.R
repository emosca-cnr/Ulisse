#' ora2enrich
#' this function translates the result of an ORA into the DOSE class "enrichResult"
#' @param ora_res ora analysis result
#' @param pvalueCutoff cut off over adjusted p value= 0.25
#' @param pAdjustMethod method for p value adjustement
#' @param qvalueCutoff cut off over q values
#' @param wb input gene list
#' @param universe universe of genes
#' @param gsl list of gene sets
#' @param organism organism
#' @param keytype keytype
#' @param ontology ontology
#' @param readable readable
#' @param min_size minimum gene set size
#' @param max_size maximum gene set size
#' @export
#' @importClassesFrom DOSE enrichResult
#' @importFrom qvalue qvalue
#' @importFrom methods new


ora2enrich <- function(ora_res=NULL, pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, wb = NULL, universe = NULL, gsl = NULL, organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN", readable = FALSE, min_size=5, max_size=500){
  
  gsl <- filter_gsl(gsl = gsl, universe = universe, min_size = min_size, max_size = max_size)
  gsl_size <- lengths(gsl)
  
  if(length(gsl)==0){
    stop("Problem with gene sets. Check identifiers and sizes.\n")
  }
  
  ora_res_gsl <- lapply(ora_res, function(x) x$id)
  gsl <- gsl[names(gsl) %in% ora_res_gsl]
  
  ans <- ora_res
  for(i in 1:length(ora_res)){
    
    geneID <- unlist(sapply(ora_res[[i]]$genes, function(x) gsub(";", "/", x)))
    
    Over <- data.frame(ID = as.character(ora_res[[i]]$id), GeneRatio = paste0(ora_res[[i]]$wbd, "/",  ora_res[[i]]$bd), BgRatio = paste0(ora_res[[i]]$wb, "/",  ora_res[[i]]$N), pvalue = ora_res[[i]]$p, p.adjust=ora_res[[i]]$p_adj, qvalue=ora_res[[i]]$q_val, Count=ora_res[[i]]$wbd, geneID=geneID, stringsAsFactors = FALSE)
    
    Over$Description <- ora_res[[i]]$description
    row.names(Over) <- as.character(Over$ID)
    
    ans[[i]] <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff,
             pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
             gene = as.character(wb[[i]]), universe = universe, geneSets = gsl,
             organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN",
             readable = FALSE)
    
    
  }
  
  return(ans)
  
  
}
