#' ulisse2enrich
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


uOra2enrich <- function(ulisse_res, pvalueCutoff = 0.25, pAdjustMethod = "BH", qvalueCutoff = 0.25, gene = NULL, universe = NULL, geneSets = NULL, organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN", readable = FALSE){

  #testing
  # ulisse_res <- read.delim("/home/bioinformatics/vnale/GR4MS/prova_GR4MS/ora/brain/bs_reactome_complete_ora_res.txt", sep=" ", stringsAsFactors = F)
  # pvalueCutoff <- 0.5
  # qvalueCutoff <- 0.25
  # pAdjustMethod <- "BH"
  #
  # load("/home/bioinformatics/vnale//GR4MS/interactomes_data/brain_graph_top1.RData")
  # universe <- igraph::V(tissue_graph)$name
  #
  # source("~/source/general/annotations/load_biosystems_collections.R")
  # bs <- load_biosystems_collections(c("bs_kegg", "bs_reactome", "bs_GO"), size = c(5, 500))
  # ms1 <- msigdbr(category = "C2", subcategory = "CP:BIOCARTA")
  # ms2 <- msigdbr(category = "H")
  # gs <- c(bs$gsl, list(biocarta=split(ms1$entrez_gene, ms1$gs_id), HALLMARKS=split(ms2$entrez_gene, ms2$gs_id)))
  # gsid2name <- rbind(unique(bs$df[, c("gsid", "db_source", "name")]),
  #                    unique(data.frame(gsid=ms1$gs_id, db_source=ms1$gs_subcat, name=ms1$gs_name, stringsAsFactors = F)),
  #                    unique(data.frame(gsid=ms2$gs_id, db_source=ms2$gs_subcat, name=ms2$gs_name, stringsAsFactors = F)))
  # rm(bs, ms1, ms2)
  # geneSets <- gs$bs_reactome
  # load("/home/bioinformatics/vnale/GR4MS/prova_GR4MS/brain/neda/top1/final_top_network.RData")
  # gene <- V(top_network)$name
  #geneID <- unlist(lapply(geneSets[match(ulisse_res$gsid, names(geneSets))], function(x) paste0(x, collapse = "/")))

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
