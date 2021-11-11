#' Plot a comparion between two ORAs
#' @param ora_res_1 result of ora_pipeline
#' @param ora_res_2 another result of ora_pipeline
#' @param p_sig thresold for significant pathways
#' @param p_max maximum p-value allowed for marginally significant pathways
#' @param dir_out output directory
#' @param ora_names project names
#' @param mar mar parameter
#' @param mgp mgp parameter
#' @param cex.axis cex.axis
#' @param use_nominal_p whether to use or not nominal p-values
#' @export

plot_ora_comparison <- function(ora_res_1=NULL, ora_res_2=NULL, p_sig=0.001, p_max=0.1, dir_out="./", ora_names=c("1", "2"), mar=c(2.5, 20, 1, 1), mgp=c(1.2, 0.3, 0), cex.axis=0.6, use_nominal_p=TRUE){
  
  p_type <- "p"
  if(!use_nominal_p){
    p_type <- "p_adj"
  }
  ora_res_merged <- merge(ora_res_1[, c("gsid", "name", "er", p_type)], ora_res_2[, c("gsid", "name", "er", p_type)], by=c("gsid", "name"), sort=F, suffixes = c(paste0("_", ora_names)))
  p_columns <- which(grepl("^p_", colnames(ora_res_merged)))
  er_columns <-  which(grepl("^er_", colnames(ora_res_merged)))
  
  ora_res_merged_filt <- ora_res_merged[apply(ora_res_merged[, p_columns], 1, function(x) any(x<p_sig)), ]
  print(nrow(ora_res_merged_filt))
  ora_res_merged_filt <- ora_res_merged_filt[order(apply(ora_res_merged_filt[, p_columns], 1, function(x) max(-log10(x)))), ]
  both_ok <- apply(ora_res_merged_filt[, p_columns], 1, function(x) all(x<p_max))
  ora_res_merged_filt$both_under_p_max <- FALSE
  ora_res_merged_filt$both_under_p_max[both_ok] <- TRUE
  
  cut_val <- c(ora_res_merged_filt[, p_columns[1]], ora_res_merged_filt[, p_columns[2]])
  cex_lev <- cut(-log10(cut_val), breaks = c(-p_max, -log10(p_max), -log10(p_sig), -log10(min(cut_val))))
  
  dir.create(dir_out, recursive = T)
  
  jpeg(paste0(dir_out, "/ora_res_dotplot.jpg"), width = 200, height = 200, units="mm", res=300)
  
  par(mar=mar)
  par(mgp=mgp)
  
  plot(ora_res_merged_filt[, er_columns[1]], 1:nrow(ora_res_merged_filt), cex=c(0.4, 1, 1.6)[as.numeric(cex_lev)[1:nrow(ora_res_merged_filt)]], pch=16, col=adjustcolor("purple", 0.8), yaxt="none", ylab="", xlab="ER", cex.axis=cex.axis)
  
  points(ora_res_merged_filt[, er_columns[2]], 1:nrow(ora_res_merged_filt), cex=c(0.4, 1, 1.6)[as.numeric(cex_lev)[(nrow(ora_res_merged_filt)+1):length(cex_lev)]], pch=16, col=adjustcolor("orange", 0.8), yaxt="none", ylab="", xlab="ER", cex.axis=cex.axis)
  
  abline(h=1:nrow(ora_res_merged_filt), col="gray", lty=3)
  
  axis(2, at=1:nrow(ora_res_merged_filt), ora_res_merged_filt$name, cex.axis=cex.axis, las=2, tick = F)
  axis(2, at=(1:nrow(ora_res_merged_filt))[both_ok], ora_res_merged_filt$name[both_ok], cex.axis=cex.axis, las=2, tick = F, font=2)
  
  legend("bottomright", pch=16, pt.cex=c(0.4,1, 1.6, 1, 1), legend =c(paste0(">", p_max), paste0("<", p_max), paste0("<", p_sig), ora_names[1], ora_names[2]), cex=cex.axis, col=c(1, 1, 1, "purple", "orange"))
  
  dev.off()
  
  return(ora_res_merged_filt)
  
}