#' Over Representation Analysis
#' @param wb list of character vectors with hits (white balls)
#' @param universe universe, character vector
#' @param wbd_min minimum number of white balls drawn
#' @param gsl named list of sets
#' @param p_adj_method p value adjustment method, see p.adjust.methods
#' @param min_size minimum gene set size
#' @param max_size maximum gene set size
#' @param out_file_prefix prefix for .xlsx and .txt output files
#' @param description optional named vector with gene set description; names must be gene set identifiers
#' @import openxlsx
#' @importFrom qvalue qvalue
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @export

ora <- function(wb=NULL, universe=NULL, gsl=NULL, p_adj_method='fdr', description=NULL, wbd_min=1, min_size=5, max_size=500, out_file_prefix="ora_res"){
  
  #gene set size
  cat("Checking gene sets\n")
  
  gsl <- filter_gsl(gsl = gsl, universe = universe, min_size = min_size, max_size = max_size)
  gsl_size <- lengths(gsl)
  
  if(length(gsl)==0){
    stop("Problem with gene sets. Check identifiers and sizes.\n")
  }
  
  if(!is.null(description)){
    description <- description[match(names(gsl), names(description))]
  }else{
    description <- setNames(names(gsl), names(gsl))
  }
  if(length(description)==0){
    description <- setNames(names(gsl), names(gsl))
  }
  
  
  cat("#{gene sets} that overlap with the universe and have size in [", min_size, ",", max_size, "]:", length(gsl), "\n")
  
  out <- wb
  
  for(i in 1:length(wb)){
    
    bb <- universe[!universe %in% wb[[i]]]
    out[[i]] <- lapply(gsl, function(x) ora1gs(wb[[i]], bb, x))
    out[[i]] <- as.data.frame(do.call(rbind, out[[i]]), stringsAsFactors = FALSE)
    
    if(wbd_min>0){
      out[[i]] <- out[[i]][out[[i]]$wbd >= wbd_min, ]
    }
    
    cat("#{gene sets} with wbd >=", wbd_min, ":", nrow(out[[i]]), "\n")
    
    out[[i]]$N <- length(wb[[i]]) + length(bb)
    out[[i]]$exp <- out[[i]]$wb * out[[i]]$bd / out[[i]]$N
    out[[i]]$id <- rownames(out[[i]])
    out[[i]]$p_adj <- p.adjust(out[[i]]$p, method = p_adj_method)
    out[[i]]$genes <-   unlist(lapply(gsl[match(out[[i]]$id, names(gsl))], function(x) paste0(sort(x[x %in% wb[[i]]]), collapse = ";")))
    
    #### FROM DOSE
    qobj <- tryCatch(qvalue(out[[i]]$p, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
    if (inherits(qobj, "qvalue")) {
      qvalues <- qobj$qvalues
    } else {
      qvalues <- NA
    }
    
    out[[i]]$q_val <- qvalues
    out[[i]]$er <- out[[i]]$wbd / out[[i]]$exp
    out[[i]]$description <- description[match(out[[i]]$id, names(description))]
    
    out[[i]] <- out[[i]][order(out[[i]]$p_adj, -out[[i]]$er), c('id', 'N', 'wb', 'bb', 'bd', 'wbd', 'exp', 'er', 'p', 'p_adj', 'q_val', 'description', 'genes')]
  }
  
  #writing output to file
  if(!is.null(out_file_prefix)){
    
    out_file_xlsx <- paste0(out_file_prefix, ".xlsx")
    cat("Writing output to", out_file_xlsx, " and ", paste0(out_file_prefix, ".[run].txt"), "...\n")
    wb <- createWorkbook()
    
    ### legend
    legend_txt <- data.frame(column=c("N", "wb", "bb", "bd", "wbd", "exp", "er", "p", "p_adj", "q_val", "genes"), description=c("all considered genes", "genes in the input set", "genes not in the input set", "genes in the pathway", "genes in the input set & in the pathway", "expected wbd in an hypergeometric experiment", "enrichment ratio", "hypergeometric p", "FDR (BH)", "FDR (qvalue)", "genes in the input set & in the pathway"), stringsAsFactors = F)
    addWorksheet(wb, "Legend")
    writeData(wb, "Legend", legend_txt)
    
    for(i in 1:length(out)){
      
      addWorksheet(wb, names(out)[i])
      writeData(wb, names(out)[i], out[[i]])
      
      write.table(out[[i]], file = paste0(out_file_prefix, ".", names(out)[i], ".txt"), row.names = F, sep="\t")
      
    }
    
    saveWorkbook(wb, out_file_xlsx, TRUE)
    
  }
  
  
  return(out)
  
}
