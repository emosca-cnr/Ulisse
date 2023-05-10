#' Function to plot functional relevance results
#' @description Graphical function to plot (relative) functional relevance results.
#' @details The function plots (relative) functional relevance results obtained with `gene_functional_relevance()` function. Each gene is plotted by interactor
#'  and functional diversity, and is colored by functional relevance score
#' @param fr the data.frame obtained as `gene_functional_relevance()` result
#' @param method = c("count", "relative") If "count" then the functional and interactor diversity and functional relevance are used, otherwise the relative ones
#' @param file_name file name of the plot, if `NULL` the plot will be returned by the function
#' @param plot_names logical, or a vector with the same length and order of the genes in `fr`. If `TRUE` all genes name are plotted,
#'  if `FALSE` no names are plotted, otherwise the passed vector is used to label genes and all of them are plotted
#' @param pal vector of three color to be used to create the gradient to color the genes according to (relative) functional relevance score. If `NULL`
#'  the Brewer "Spectral" palette is used
#' @param jitter_width,jitter_height values to control point jitter in the plot. See ggplot2::geom_jitter() for details
#' @param width,height,res,units values used to save the plot in jpeg format, see grDevices::jpeg() for details
#' @return If `file_name` The function returns the ggplot object
#' @import plotrix
#' @import ggplot2
#' @import ggrepel
#' @importFrom grDevices dev.off jpeg
#' @importFrom scales rescale
#' @export



plot_functional_relevance <- function(fr, method = "count", file_name = NULL, plot_names =T, pal = NULL, jitter_width = 0.07, jitter_height = 0.07,
                                      width = 200, height = 200,
                                      res = 300, units = "mm") {
  
  tab.out <- fr
  if(method == "count") {
    tab.out <- tab.out[, c("gene", "functional_diversity", "interactor_diversity", "functional_relevance")]
  } else if(method == "relative") {
    tab.out <- tab.out[, c("gene", "relative_functional_diversity", "relative_interactor_diversity", "relative_functional_relevance")]
    colnames(tab.out) <- c("gene", "functional_diversity", "interactor_diversity", "functional_relevance")
  }
  
  if(is.logical(plot_names)) {
    if(plot_names) {
      tab.out$label <- tab.out$gene
    } else {
      tab.out$label <- ""
    }
  } else {
    tab.out$label <- plot_names
  }
  
  Mf <- max(tab.out$functional_relevance)
  Mf <- signif(Mf, digits = 3)
  mf <- min(tab.out$functional_relevance)
  mf <- signif(mf, digits = 3)
  if(is.null(pal)) {
    pal <- RColorBrewer::brewer.pal(n = 5,  name = "Spectral")[c(5,3,1)]
  }
  
  
  p <- ggplot(tab.out, aes_string(y="functional_diversity", x="interactor_diversity", label = "label")) + 
    geom_jitter(aes_string(fill = "functional_relevance"), width = 0.07, height = 0.07, shape = 21, color = "gray75",stroke = 0.2) +
    geom_text_repel(fontface="italic", max.overlaps = 2000)+
    theme_light() +
    labs(fill = expression(italic("f"))) +
    geom_abline(intercept = 0, color = "gray65", linetype = "dashed") +
    scale_fill_gradientn(colours = pal,
                         values = scales::rescale(c(mf,0,Mf)), breaks = c(mf,0,Mf), limits = c(mf-0.1, Mf + 0.1)) +
    ylab(bquote(n[f])) + xlab(bquote(n[d]))
  
  if(!is.null(file_name)) {
    jpeg(file_name, width = width, height = height,
         res = res, units = units)
  }
  
  print(p)
  
  if(!is.null(file_name)) {
    dev.off()
  }
  
  if(is.null(file_name)) {
    return(p)
  }
}
