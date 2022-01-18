#' Obtaining data from TCGA for vignette
#' @param project name of the TCGA project
#' @param data.category category of the data
#' @param data.type type of data
#' @param sample.type type of sample
#' @param workflow.type type of workflow
#' @param summarizedExperiment logical indicating if the data should be imported as a summarized experiment
#' @param save logical, save or not the downloaded data locally
#' @import TCGAbiolinks
#' @export

get_TCGA <- function(project, data.category, data.type, sample.type, workflow.type, summarizedExperiment, save) {
  query <- TCGAbiolinks::GDCquery(project = project, data.category = data.category,
                    data.type = data.type, sample.type = sample.type,
                    workflow.type = workflow.type
                    )
  TCGAbiolinks::GDCdownload(query = query) #have to be done only the first time
  data <- TCGAbiolinks::GDCprepare(query, summarizedExperiment = summarizedExperiment, save = save)
  return(data)

}