#' Gene Sert Enrichment Analysis
#' @param rl numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory
#' @param gsl named list of gene sets
#' @param k integer, number of permutations
#' @param min_size minimum gene set size
#' @param max_size maximum gene set size
#' @param decreasing TRUE/FALSE vector that specifies whether to order each column of rl decreasingly or not; must be of length equal to `ncol(rl)`. If NULL, all columns will be ranked in decreasing order
#' @param mc_cores_path number of cores to use for parallel calculation of gene set lists; the total number of cpu used will be mc_cores_path x mc_cores_perm
#' @param mc_cores_perm number of cores to use for parallel calculation of ranked list permutations; the total number of cpu used will be mc_cores_path x mc_cores_perm
#' @param description optional named vector with gene set description; names must be gene seet identifiers
#' @param out_file_prefix prefix for .xlsx and .txt output files
#' @param min.k minimum number of permutations to obtain valid permutation-based statistics
#' @param min_tags minimum number of tags to consider the ES; gene sets with tags < min_tags will be excluded
#' @import parallel openxlsx
#' @importFrom qvalue qvalue
#' @importFrom utils write.table
#' @return data.frame with: es, enrichment score; nes normalized enrichment score; nperm, number of permutations actually used; p-value, empirical p-value; adjusted p-value, BH FDR; q_val: q-value estimnated from p-values using qvalue package; FDR q-value, empirical FDR; tags, leading edge size; tags_perc, leading edge size percent over gene set; list_top, rank of the ES; list_top_perc, rank of the ES percent over full ranked list; lead_edge, signal strength; lead_edge_subset, gene names of the leading edge
#' @export

gsea <- function(rl=NULL, gsl=NULL, k=99, min_size=5, max_size=500, min_tags=3, decreasing=NULL, mc_cores_path=1, mc_cores_perm=1, description=NULL, out_file_prefix="gsea_res", min.k=20){
  
  #checks
  if(!is.matrix(rl) | !is.numeric(rl)){
    stop("rl must be a numeric matrix")
  }
  
  if(is.null(decreasing)){
    decreasing <- rep(TRUE, ncol(rl))
  }
  
  if(length(decreasing) != ncol(rl)){
    stop("length(decreasing) must be equal to ncol(rl)")
  }
  
  #gene set size
  cat("Checking gene sets\n")
  
  gsl <- filter_gsl(gsl = gsl, universe =  rownames(rl), min_size = min_size, max_size = max_size)
  gsl_size <- lengths(gsl)
  
  if(length(gsl)==0){
    stop("Problem with gene sets. Check identifiers and sizes.\n")
  }
  
  cat("#{gene sets} with genes occurring in the ranked lists and size in [", min_size, ",", max_size, "]:", length(gsl), "\n")
  
  
  if(!is.null(description)){
    description <- description[match(names(gsl), names(description))]
  }else{
    description <- setNames(names(gsl), names(gsl))
  }
  if(length(description)==0){
    description <- setNames(names(gsl), names(gsl))
  }
  
  #create the list of ranked vectors
  cat("Decreasing:", decreasing, "\n")
  rll <- vector('list', ncol(rl))
  names(rll) <- colnames(rl)
  for(i in 1:length(rll)){
    #rll[[i]] <- sort(array(rl[, i], dimnames = list(rownames(rl))), decreasing = decreasing[i])
    rll[[i]] <- sort(setNames(rl[, i], rownames(rl)), decreasing = decreasing[i])
    if(length(unique(rll[[i]])) != length(rll[[i]])){
      cat("Found ties in ranked list ", names(rll)[i], "!!!\n")
    }
  }
  
  #permutation of gene ids
  cat('generating', k, 'permutations\n')
  x_perm <- lapply(1:k, function(x) sample(rownames(rl), nrow(rl)))
  
  #real es
  cat("ES of input data...\n")
  real_es_data <- lapply(gsl, function(x) lapply(rll, function(y) es(which(names(y) %in% x), y)))
  real_es <- do.call(rbind, lapply(real_es_data, function(x) unlist(lapply(x, function(y) y$es))))
  
  leading_edge <- vector("list", length(rll))
  names(leading_edge) <- names(rll)
  for(i in 1:length(leading_edge)){
    leading_edge[[i]] <- do.call(rbind, lapply(real_es_data, function(x) x[[i]][, -1])) #rm the enrichment score
  }
  
  #permutations
  cat("ES of permutations...\n")
  res <- mclapply(gsl, function(x) do.call(rbind, mclapply(x_perm, function(y) calc_gs_perm(rll, y, x), mc.cores=mc_cores_perm)), mc.cores = mc_cores_path)
  
  temp <- vector('list', length(rll))
  names(temp) <- colnames(rl)
  for(i in 1:length(rll)){
    temp[[i]] <- cbind(real_es[, i], do.call(rbind, lapply(res, function(x) x[, i])))
  }
  res <- temp
  rm(temp)
  
  
  #statistics
  print("calculating statistics...")
  
  #the first column is the real value, and it is included in the calculation of p
  #if ES* > 0 -> p = # (ESp >= ES*) / (positive null)
  #if ES* < 0 -> p = # (ESp <= ES*) / (negative null)
  out <- res
  
  for(i in 1:length(rll)){
    
    cat("Proceeding with", names(rll)[i], "...\n")
    n_pos_perm <- rowSums(res[[i]]>0)
    n_neg_perm <- rowSums(res[[i]]<0)
    
    p_val <- apply(res[[i]], 1, function(x) ifelse(x[1] >= 0, sum(x >= x[1]) / length(x[x>=0]), sum(x <= x[1]) / length(x[x<=0])))
    
    idx_na <- which(res[[i]][, 1] == 0)
    if(length(idx_na)>0){
      cat("\tES==0:\n")
      cat("\t", rownames(res[[i]])[idx_na], "\n")
      p_val[idx_na] <- NA ### 
    }
    
    idx_na <- which(res[[i]][, 1] > 0 & n_pos_perm < min.k)
    if(length(idx_na)>0){
      cat("\tES>0 but less than", min.k, " positive ES in permutations\n")
      cat("\t", rownames(res[[i]])[idx_na], "\n")
      p_val[idx_na] <- NA ### 
    }
    
    idx_na <- which(res[[i]][, 1] < 0 & n_neg_perm < min.k)
    if(length(idx_na)>0){
      cat("\tES<0 but less than", min.k, " negative ES in permutations\n")
      cat("\t", rownames(res[[i]])[idx_na], "\n")
      p_val[idx_na] <- NA ### 
    }
    
    idx_na <- which(leading_edge[[i]]$tags < min_tags)
    if(length(idx_na)>0){
      cat("\tES is based on less than", min_tags, "\n")
      cat("\t", rownames(res[[i]])[idx_na], "\n")
      res[[i]][idx_na, 1] <- NA ### ES <- NA
      p_val[idx_na] <- NA ### 
      
      if(all(is.na(p_val))){
        cat("\t All ES are based on less than", min_tags, "\n")
        return(NULL)
      }
      
      #remove ES with tags < min_tags
      res[[i]] <- res[[i]][-idx_na, ]
      p_val <- p_val[-idx_na]
      n_pos_perm <- n_pos_perm[-idx_na]
      n_neg_perm <- n_neg_perm[-idx_na]
      
    }
    
    #normalized ES
    means <- t(apply(res[[i]], 1, function(x) c(mean(x[x>0]), abs(mean(x[x<0]))))) #positive, negative
    means[is.nan(means)] <- NA #NaN values are caused by the absence of any positive or negative value
    nes <- res[[i]] / means[, 1]
    nes_neg <- res[[i]] / means[, 2]
    nes[res[[i]] < 0] <- nes_neg[res[[i]] < 0]
    nes[is.nan(nes)] <- NA #NaN values are caused by 0/0
    rm(means, nes_neg)
    
    nes[nes > 0 & n_pos_perm < min.k] <- NA
    nes[nes < 0 & n_neg_perm < min.k] <- NA
    
    #calculate FDR
    all_nes <- as.numeric(nes)
    nes_pos <- all_nes[which(all_nes > 0)]
    nes_neg <- all_nes[which(all_nes < 0)]
    n_nes_pos <- length(nes_pos)
    n_nes_neg <- length(nes_neg)
    
    real_nes <- nes[, 1]
    real_nes <- real_nes[!is.na(real_nes)]
    real_nes_pos <- real_nes[which(real_nes>0)]
    real_nes_neg <- real_nes[which(real_nes<0)]
    n_real_nes_pos <- length(real_nes_pos)
    n_real_nes_neg <- length(real_nes_neg)
    
    #FDR: NES* > 0: fdrq = [#(all positive NESp >= NES*) / #(all positive NESp)] / [ #(all positive NES* >= NES*) / (all positive NES*) ]
    #FDR: NES* < 0: fdrq = [#(all negative NESp <= NES*) / #(all negative NESp)] / [ #(all negative NES* <= NES*) / (all negative NES*) ]
    #fdrq <- sapply(nes[, 1], function(x) ifelse(x>0, sum(all_nes >= x, na.rm = T) / n_nes_pos, sum(all_nes <= x, na.rm = T) / n_nes_neg))
    #fdrq <- fdrq / sapply(nes[, 1], function(x) ifelse(x>0, sum(nes[, 1] >= x, na.rm = T) / n_real_nes_pos, sum(nes[, 1] <= x, na.rm = T) / n_real_nes_neg))
    
    fdrq <- sapply(nes[, 1], function(x) ifelse(x>0, sum(nes_pos >= x, na.rm = T) / n_nes_pos, sum(nes_neg <= x, na.rm = T) / n_nes_neg))
    fdrq <- fdrq / sapply(nes[, 1], function(x) ifelse(x>0, sum(real_nes_pos >= x, na.rm = T) / n_real_nes_pos, sum(real_nes_neg <= x, na.rm = T) / n_real_nes_neg))
    
    rm(all_nes)
    
    fdrq[fdrq>1] <- 1
    idx_replace <- which(fdrq<p_val)
    fdrq[idx_replace] <- p_val[idx_replace]
    
    
    #### FROM DOSE
    qobj <- tryCatch(qvalue(p_val, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
    if (inherits(qobj, "qvalue")) {
      qvalues <- qobj$qvalues
    } else {
      qvalues <- NA
    }
    
    nperm <- n_pos_perm
    nperm[which(res[[i]][, 1]<0)] <- n_neg_perm[which(res[[i]][, 1]<0)]
    nperm[which(res[[i]][, 1] == 0)] <- NA
    
    #out table
    out[[i]] <- data.frame(id=rownames(res[[i]]), size=gsl_size[match(rownames(res[[i]]), names(gsl_size))], es=res[[i]][, 1], nes=nes[, 1], nperm=nperm, p_val=p_val, adj_p_val=stats::p.adjust(p_val, method='fdr'), q_val=qvalues, FDRq=fdrq, description=description[match(rownames(res[[i]]), names(description))], stringsAsFactors=FALSE)
    
    out[[i]] <- merge(out[[i]], leading_edge[[i]], by.x=1, by.y=0, sort=F)
    
    ##order by significance
    out[[i]] <- out[[i]][order(out[[i]]$FDRq, -abs(out[[i]]$nes)), ]
    
  }
  
  #writing output to file
  if(!is.null(out_file_prefix)){
    
    out_file_xlsx <- paste0(out_file_prefix, ".xlsx")
    cat("Writing output to", out_file_xlsx, " and ", paste0(out_file_prefix, ".[run].txt"), "...\n")
    
    wb <- createWorkbook()
    
    legend_txt <- data.frame(column=c("size", "es", "nes", "nperm", "p_val", "adj_p_val", "q_val", "FDRq", "tags", 	"tags_perc", "list_top", "list_top_perc", "lead_edge", "lead_edge_subset"), description=c("gene set size", "enrichment score", "normalized enriched score", "number of permutations actually used", "empirical p-value", "FDR (BH)", "FDR (qvalue)", "FDR (empirical)", "leading edge size", "leading edge size percent over gene set", "rank of the ES", "rank of the ES percent over full ranked list", "signal strength", "gene names of the leading edge"), stringsAsFactors = F)
    addWorksheet(wb, "Legend")
    writeData(wb, "Legend", legend_txt)
    
    for(i in 1:length(out)){
      
      addWorksheet(wb, names(out)[i])
      writeData(wb, names(out)[i], out[[i]])
      
      write.table(out[[i]], file = paste0(out_file_prefix, ".", names(out)[i], ".txt"), row.names = F, sep="\t")
      
    }
    
    saveWorkbook(wb, out_file_xlsx, TRUE)
  }
  
  cat("done\n")
  
  return(out)
  
}
