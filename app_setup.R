library(shiny)
library(plotly)

GSE121212_data <- readRDS("./GSE121212/GSE121212_object.rds")
condition_AD_vs_CTRL <- readRDS( "./GSE121212/resDEA_AD_VS_CTRL.rds")
genelist <- GSE121212_data$expression_data$gene 



makeLongExpData <- function(inData, inSelectedGenes, inMetadata ) {
  inDataFiltered <- inData |>
    dplyr::filter(gene %in% inSelectedGenes)
  
  longExprData <- inDataFiltered |>
    reshape2::melt(id.vars = "gene",
                   variable.name = "sample",
                   value.name = "expression")
  
  longExprData <- longExprData |>
    dplyr::left_join(inMetadata[c("title", "skin type:ch1", "patient's condition:ch1" )], by = c("sample"="title"))
  
  return(longExprData)
  
}
