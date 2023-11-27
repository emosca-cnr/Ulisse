#' Create a gseaResult instance from Ulisse GSEA result
#' @param gsea_res output of function gsea()
#' @param rl the input used to obtain gsea_res; numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory;
#' @param gsl named list of gene sets used to obtain gsea_res
#' @param min_size minimum gene set size
#' @param max_size maximum gene set size
#' @param decreasing TRUE/FALSE vector that specifies whether to order each column of rl decreasingly or not; must be of length equal to `ncol(rl)`. If NULL, all #' @export
#' @importClassesFrom DOSE gseaResult
#' @importFrom qvalue qvalue
#' @importFrom methods new

gsea2enrich <- function(gsea_res=NULL, rl=NULL, gsl = NULL, min_size=5, max_size=500, decreasing=NULL){
  
  #from https://rdrr.io/github/GuangchuangYu/DOSE/src/R/gsea.R
  gseaScores <- function(geneList=NULL, geneSet=NULL, exponent=1, fortify=FALSE) {
  
    geneSet <- intersect(geneSet, names(geneList))
    
    N <- length(geneList)
    Nh <- length(geneSet)
    
    Phit <- Pmiss <- numeric(N)
    hits <- names(geneList) %in% geneSet ## logical
    
    Phit[hits] <- abs(geneList[hits])^exponent
    NR <- sum(Phit)
    Phit <- cumsum(Phit/NR)
    
    Pmiss[!hits] <-  1/(N-Nh)
    Pmiss <- cumsum(Pmiss)
    
    runningES <- Phit - Pmiss
    
    ## ES is the maximum deviation from zero of Phit-Pmiss
    max.ES <- max(runningES)
    min.ES <- min(runningES)
    if( abs(max.ES) > abs(min.ES) ) {
      ES <- max.ES
    } else {
      ES <- min.ES
    }
    
    df <- data.frame(x=seq_along(runningES),
                     runningScore=runningES,
                     position=as.integer(hits)
    )
    
    if(fortify==TRUE) {
      return(df)
    }
    
    df$gene = names(geneList)
    res <- list(ES=ES, runningES = df)
    
    return(res)
  }
  
  leading_edge <- function(observed_info) {
    core_enrichment <- lapply(observed_info, function(x) {
      runningES <- x$runningES
      runningES <- runningES[runningES$position == 1,]
      ES <- x$ES
      if (ES >= 0) {
        i <- which.max(runningES$runningScore)
        leading_gene <- runningES$gene[1:i]
      } else {
        i <- which.min(runningES$runningScore)
        leading_gene <- runningES$gene[-c(1:(i-1))]
      }
      return(leading_gene)
    })
    
    rank <- sapply(observed_info, function(x) {
      runningES <- x$runningES
      ES <- x$ES
      if (ES >= 0) {
        rr <- which.max(runningES$runningScore)
      } else {
        i <- which.min(runningES$runningScore)
        rr <- nrow(runningES) - i + 1
      }
      return(rr)
    })
    
    tags <- sapply(observed_info, function(x) {
      runningES <- x$runningES
      runningES <- runningES[runningES$position == 1,]
      ES <- x$ES
      if (ES >= 0) {
        i <- which.max(runningES$runningScore)
        res <- i/nrow(runningES)
      } else {
        i <- which.min(runningES$runningScore)
        res <- (nrow(runningES) - i + 1)/nrow(runningES)
      }
      return(res)
    })
    
    ll <- sapply(observed_info, function(x) {
      runningES <- x$runningES
      ES <- x$ES
      if (ES >= 0) {
        i <- which.max(runningES$runningScore)
        res <- i/nrow(runningES)
      } else {
        i <- which.min(runningES$runningScore)
        res <- (nrow(runningES) - i + 1)/nrow(runningES)
      }
      return(res)
    })
    
    N <- nrow(observed_info[[1]]$runningES)
    setSize <- sapply(observed_info, function(x) sum(x$runningES$position))
    signal <- tags * (1-ll) * (N / (N - setSize))
    
    tags <- paste0(round(tags * 100), "%")
    ll <- paste0(round(ll * 100), "%")
    signal <- paste0(round(signal * 100), "%")
    leading_edge <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)
    
    res <- list(rank = rank,
                tags = tags,
                list = ll,
                signal = signal,
                leading_edge = leading_edge,
                core_enrichment = core_enrichment)
    return(res)
  }  

  params <- list(pvalueCutoff = 1,
                 nPerm = NA,
                 pAdjustMethod = "BH",
                 exponent = 1,
                 minGSSize = NA,
                 maxGSSize = NA
  )
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
  gsl <- filter_gsl(gsl = gsl, universe =  rownames(rl), min_size = min_size, max_size = max_size)
  gsl_size <- lengths(gsl)

  #create the list of ranked vectors
  cat("Decreasing:", decreasing, "\n")
  rll <- vector('list', ncol(rl))
  names(rll) <- colnames(rl)
  ans <- rll
  
  for(i in 1:length(rll)){
    
    rll[[i]] <- sort(setNames(rl[, i], rownames(rl)), decreasing = decreasing[i])
    if(length(unique(rll[[i]])) != length(rll[[i]])){
      cat("Found ties in ranked list ", names(rll)[i], "!!!\n")
    }
    
    ans[[i]] <- data.frame(
      ID = as.character(gsea_res[[i]]$id),
      Description = as.character(gsea_res[[i]]$description),
      setSize = gsl_size[match(gsea_res[[i]]$id, names(gsl_size))],
      enrichmentScore = gsea_res[[i]]$es,
      NES = gsea_res[[i]]$nes,
      pvalue = gsea_res[[i]]$p_val,
      p.adjust = gsea_res[[i]]$adj_p_val,
      qvalue = gsea_res[[i]]$q_val,
      stringsAsFactors = FALSE
    )
    
    idx <- order(ans[[i]]$p.adjust, -abs(ans[[i]]$NES), decreasing = FALSE)
    ans[[i]] <- ans[[i]][idx, ]
    
    rownames(ans[[i]]) <- ans[[i]]$ID
    
    observed_info <- lapply(gsl[ans[[i]]$ID], function(gs) gseaScores(geneSet=gs, geneList=rll[[i]]))
    
    ledge <- leading_edge(observed_info)
    
    ans[[i]]$rank <- ledge$rank
    ans[[i]]$leading_edge <- ledge$leading_edge
    ans[[i]]$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')
    
    ans[[i]] <- new("gseaResult",
        result     = ans[[i]],
        geneSets   = gsl,
        geneList   = rll[[i]],
        permScores = matrix(),
        params     = params,
        readable   = FALSE
    )
    
  }
  
  return(ans)
  
}