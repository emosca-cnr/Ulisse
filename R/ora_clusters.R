#' Over Representation Analysis on all clusters of single cell data
#' @description This function can be used to perform ORA test on multiple clusters/cell-types or conditions
#' @param wb named list with hits (white balls) vector for each cluster/condition
#' @param universe named list with universe genes for each cluster/condition. Will then be used to obtain black balls ("other elements")
#' @param gsl named list of sets
#' @param wb_min valued used to filter gsl used for testing. If equal to 0, all sets will be tested, otherwise only the ones that 
#' contains at least `wb_min` `wb` genes will be tested. Useful to limit number of tested sets and avoid testing sets that do not contains `wb` genes
#' @param p_adj_method p value adjustment method, see p.adjust.methods
#' @param mc_cores number of cores used to parallelize ORA tests
#' @return The function return a named list with a result table for each cluster/condition. In each table the columns are
#' \itemize{
#'  \item id: name of the set tested
#'  \item N: total number of gene tested
#'  \item wb: number of hit genes
#'  \item bb: number of other elements present in the test (`N` - `wb`)
#'  \item bd: number of genes in the set tested
#'  \item wbd: number of hit genes present in the set tested
#'  \item exp:
#'  \item er:
#'  \item p: p-value
#'  \item p_adj: corrected p-value calculated by using `p_adj_method`
#'  \item q_val: q-value calculated using bootstrap method (see qvalue::qvalue for further details)
#'  \item gene: list of `wbd` genes separated by `;`
#' }
#' @importFrom qvalue qvalue
#' @import parallel
#' @export

ora_clusters <- function(wb, universe, gsl, wb_min = 1,  p_adj_method='fdr', mc_cores = 2){
  bb <- lapply(1:length(universe), function(x) {
    tmp <- universe[[x]]
    tmp <- tmp[!tmp %in% wb[[x]]]
    return(tmp)
  })
  names(bb) <- names(universe)
  
  gsl_all <- lapply(1:length(wb), function(x) {
    tmp <- lengths(lapply(gsl, function(j) j <- j[j %in% wb[[x]]]))
    gsl_out <- gsl[tmp >= wb_min]
    return(gsl_out)
  })
  
  names(gsl_all) <- names(wb)
  
  res <- mclapply(names(wb), function(x) {
    wb_tmp <- wb[[x]]
    bb_tmp <- bb[[x]]
    gsl_tmp <- gsl_all[[x]]
    out <- lapply(gsl_tmp, function(x) ora1gs(wb_tmp, bb_tmp, x))
    
    out <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
    out$N <- length(wb_tmp) + length(bb_tmp)
    out$exp <- out$wb * out$bd / out$N
    out$id <- rownames(out)
    out$p_adj <- stats::p.adjust(out$p, method = p_adj_method)
    out$q_val <- qvalue::qvalue(p=out$p, lambda=0.05, pi0.method="bootstrap")$qvalues
    out$er <- out$wbd / out$exp
    tmp <- lapply(gsl_tmp, function(x) x <- paste(sort(x[which(x %in% wb_tmp)]), collapse = ";"))
    tmp <- unlist(tmp)
    out$gene <- tmp[out$id]
    
    return(out[, c('id', 'N', 'wb', 'bb', 'bd', 'wbd', 'exp', 'er', 'p', 'p_adj', 'q_val', 'gene')])
  }, mc.cores = mc_cores)
  names(res) <- names(wb)
  return(res)
  
}
