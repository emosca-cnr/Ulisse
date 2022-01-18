#' Plot a comparion between two ORAs
#' @param ora_res_list list of results of ora_pipeline
#' @param p_sig thresold for significant pathways
#' @param p_max maximum p-value allowed for marginally significant pathways
#' @param dir_out output directory
#' @param mar mar parameter
#' @param mgp mgp parameter
#' @param cex.axis cex.axis
#' @param use_nominal_p whether to use or not nominal p-values
#' @param col_pal color palette function
#' @return data.frame with merged results only for pathways that satisfy p_sig in at least one condition
#' @importFrom viridis turbo viridis
#' @export

plot_ora_comparison <- function(ora_res_list=NULL, p_sig=0.001, p_max=0.1, dir_out="./", mar=c(2.5, 20, 1, 1), mgp=c(1.2, 0.3, 0), cex.axis=0.6, use_nominal_p=TRUE, col_pal=NULL){

  p_type <- "p"
  if(!use_nominal_p){
    p_type <- "p_adj"
  }
  if(is.null(names(ora_res_list))){
    ora_names <- paste0("id", 1:length(ora_res_list))
  }else{
    ora_names <- names(ora_res_list)
  }


  ora_res_merged <- merge(ora_res_list[[1]][, c("gsid", "name", "er", p_type)], ora_res_list[[2]][, c("gsid", "name", "er", p_type)], by=c("gsid", "name"), sort=F, suffixes = c(paste0("_", ora_names[1:2])))
  if(length(ora_res_list) > 2){
    for(i in 3:length(ora_res_list)){
      ora_res_merged <- merge(ora_res_merged, ora_res_list[[i]][, c("gsid", "name", "er", p_type)], by=c("gsid", "name"), sort=F)
    }
  }
  colnames(ora_res_merged)[match(c("er", "p"), colnames(ora_res_merged))] <- c(paste0(c("er", "p"), "_", ora_names[i]))

  p_columns <- which(grepl("^p_", colnames(ora_res_merged)))
  er_columns <-  which(grepl("^er_", colnames(ora_res_merged)))

  ora_res_merged_filt <- ora_res_merged[apply(ora_res_merged[, p_columns], 1, function(x) any(x<p_sig)), ]
  print(nrow(ora_res_merged_filt))
  ora_res_merged_filt <- ora_res_merged_filt[order(apply(ora_res_merged_filt[, p_columns], 1, function(x) max(-log10(x)))), ]
  all_ok <- apply(ora_res_merged_filt[, p_columns], 1, function(x) all(x<p_max))
  ora_res_merged_filt$all_under_p_max <- FALSE
  ora_res_merged_filt$all_under_p_max[all_ok] <- TRUE

  p_val <- unlist(ora_res_merged_filt[, p_columns])
  p_fact <- cut(-log10(p_val), breaks = c(-p_max, -log10(p_max), -log10(p_sig), -log10(min(p_val)))) #3 values

  er_val <- unlist(ora_res_merged_filt[, er_columns])
  er_fact <- cut(er_val, breaks = 4) #4 values

  dir.create(dir_out, recursive = T)

  if(is.null(col_pal)){
    col_pal1 <- viridis::turbo(length(ora_res_list)) #conditions
    col_pal2 <- viridis::viridis(3) #p_values
  }

  grDevices::jpeg(paste0(dir_out, "/ora_res_dotplot1.jpg"), width = 200, height = 200, units="mm", res=300)

  graphics::par(mar=mar)
  graphics::par(mgp=mgp)
  graphics::layout(matrix(1:2, ncol = 2), widths = c(0.9, 0.1))

  plot(ora_res_merged_filt[, er_columns[1]], 1:nrow(ora_res_merged_filt), cex=c(0.4, 1, 1.6)[as.numeric(p_fact)[1:nrow(ora_res_merged_filt)]], pch=16, col=grDevices::adjustcolor(col_pal1[1], 0.8), yaxt="none", ylab="", xlab="ER", cex.axis=cex.axis, xlim = c(min(ora_res_merged_filt[, er_columns]), max(ora_res_merged_filt[, er_columns])))

  for(i in 2:length(ora_res_list)){
    graphics::points(ora_res_merged_filt[, er_columns[i]], 1:nrow(ora_res_merged_filt), cex=c(0.4, 1, 1.6)[as.numeric(p_fact)[(((i-1)*nrow(ora_res_merged_filt))+1):(nrow(ora_res_merged_filt)*i)]], pch=16, col=grDevices::adjustcolor(col_pal1[i], 0.8), yaxt="none", ylab="", xlab="ER", cex.axis=cex.axis)
  }

  graphics::abline(h=1:nrow(ora_res_merged_filt), col="gray", lty=3)

  graphics::axis(2, at=1:nrow(ora_res_merged_filt), ora_res_merged_filt$name, cex.axis=cex.axis, las=2, tick = F)
  graphics::axis(2, at=(1:nrow(ora_res_merged_filt))[all_ok], ora_res_merged_filt$name[all_ok], cex.axis=cex.axis, las=2, tick = F, font=2)

  graphics::par(mar=c(0, 0, 0, 1))
  graphics::plot.new()
  graphics::legend("center", pch=16, pt.cex=c(0.4, 1, 1.6, rep(1, length(ora_res_list))), legend =c(paste0(">", p_max), paste0("<", p_max), paste0("<", p_sig), ora_names), cex=cex.axis, col=c(1, 1, 1, col_pal1))

  grDevices::dev.off()


  grDevices::jpeg(paste0(dir_out, "/ora_res_dotplot2.jpg"), width = 200, height = 200, units="mm", res=300)

  graphics::par(mar=mar)
  graphics::par(mgp=mgp)
  graphics::layout(matrix(1:2, ncol = 2), widths = c(0.85, 0.15))

  plot(rep(1, nrow(ora_res_merged_filt)), 1:nrow(ora_res_merged_filt), cex=c(1, 1.3, 1.6, 2)[as.numeric(er_fact)[1:nrow(ora_res_merged_filt)]], pch=16, col=grDevices::adjustcolor(col_pal2[as.numeric(p_fact)[1:nrow(ora_res_merged_filt)]], 0.8), xaxt="none", yaxt="none", ylab="", xlab="", cex.axis=cex.axis, xlim = c(0.5, length(ora_res_list)+0.5))

  for(i in 2:length(ora_res_list)){
    graphics::points(rep(i, nrow(ora_res_merged_filt)), 1:nrow(ora_res_merged_filt), cex=c(1, 1.3, 1.6, 2)[as.numeric(er_fact)[(((i-1)*nrow(ora_res_merged_filt))+1):(nrow(ora_res_merged_filt)*i)]], pch=16, col=grDevices::adjustcolor(col_pal2[as.numeric(p_fact)[(((i-1)*nrow(ora_res_merged_filt))+1):(nrow(ora_res_merged_filt)*i)]], 0.8), cex.axis=cex.axis)
  }

  graphics::abline(h=1:nrow(ora_res_merged_filt), col="gray", lty=3)
  graphics::abline(v=1:length(ora_res_list), col="gray", lty=3)

  graphics::axis(1, at=1:length(ora_res_list), ora_names, cex.axis=cex.axis, las=2, tick = F)
  graphics::axis(2, at=1:nrow(ora_res_merged_filt), ora_res_merged_filt$name, cex.axis=cex.axis, las=2, tick = F)
  graphics::axis(2, at=(1:nrow(ora_res_merged_filt))[all_ok], ora_res_merged_filt$name[all_ok], cex.axis=cex.axis, las=2, tick = F, font=2)

  graphics::par(mar=c(mar[1], 0, mar[3], 1))
  graphics::plot.new()
  graphics::legend("center", pt.cex=c(1, 1.3, 1.6, 2), legend = levels(er_fact), cex=cex.axis, pch=1)
  graphics::legend("bottom", pch=16, legend=c(paste0(">", p_max), paste0("<", p_max), paste0("<", p_sig)), cex=cex.axis, col=col_pal2)

  grDevices::dev.off()


  return(ora_res_merged_filt)

}
