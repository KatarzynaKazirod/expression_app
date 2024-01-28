source("./app_setup.R")

ui <- fluidPage(
  titlePanel("GSE121212 expression data"),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        'geneselector',
        'Select genes',
        choices = NULL,
        selected = NULL,
        multiple = TRUE
        
      ),
      br(),
      numericInput(
        "fdr_cutoff",
        "Corrected p-value cutoff:",
        value = 0.05,
        min = min(condition_AD_vs_CTRL$padj),
        max = max(condition_AD_vs_CTRL$padj),
        step = 0.01
      ),
      br(),
      numericInput(
        "lfc_cutoff",
        "Absolute value of Log Fold Change cutoff:",
        value = 1,
        min = 0,
        max = max(condition_AD_vs_CTRL$log2FoldChange),
        step = 0.5
      )
      
    ),
    
    
    mainPanel(br(),
              tabsetPanel(
                type = "tabs",
                
                tabPanel("Raw expression", plotlyOutput("plotRawExpr")),
                tabPanel("Atopic Dermatitis vs Healthy", plotlyOutput("plotVolcano"))
                
              ))
  ),
  style = "padding: 50px;"
)

server <- function(input, output, session) {

    updateSelectizeInput(session, "geneselector", selected = character(0), choices = genelist, server = TRUE)
  
  longExprData <-  reactive({
    selected_genes <- input$geneselector
    genesExpr <-
      makeLongExpData(
        inData = GSE121212_data$expression_data,
        inSelectedGenes = selected_genes,
        inMetadata = GSE121212_data$metadata
      )
    return(genesExpr)
    
  })
  

    boxplotExpr <- reactive({
      longExprData <- longExprData()
      
      boxplot <-
        plotly::plot_ly(
          data = longExprData,
          x = ~ gene,
          y = ~ expression,
          type = "box"
        )
      
    })  
      
  output$plotRawExpr <- renderPlotly({
    boxplotExpr()

  })

  
  
  volcanoPlot <- reactive({
    selected_genes <- input$geneselector
    fdr_cutoff <- input$fdr_cutoff
    lfc_cutoff <- input$lfc_cutoff
    
    
    volcanoPlot <-
      plotly::plot_ly(
        condition_AD_vs_CTRL,
        x = ~ log2FoldChange,
        y = ~ -log10(padj),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = ifelse(condition_AD_vs_CTRL$gene %in% selected_genes, "red", "grey"),
          opacity = 0.6
        ),
        text = ~ gene
      )

    
    volcanoPlot <- volcanoPlot |>
      add_lines(
        y = -log10(fdr_cutoff),
        x = range(condition_AD_vs_CTRL$log2FoldChange),
        line = list(color = "darkgreen",
                    dash = "dot"),
        inherit = FALSE,
        showlegend = FALSE
      ) |>
      add_lines(
        y = range(-log10(condition_AD_vs_CTRL$padj)),
        x = lfc_cutoff,
        line = list(color = "darkgreen",
                    dash = "dot"),
        inherit = FALSE,
        showlegend = FALSE
      ) |>
      add_lines(
        y = range(-log10(condition_AD_vs_CTRL$padj)),
        x = -lfc_cutoff,
        line = list(color = "darkgreen",
                    dash = "dot"),
        inherit = FALSE,
        showlegend = FALSE
      )
    volcanoPlot
    
  }) 
  

  
output$plotVolcano <- renderPlotly({

  volcanoPlot <- volcanoPlot()
  
  
  volcanoPlot <- volcanoPlot |>
    plotly::layout(
      title = "",
      xaxis = list(title = "log2 Fold Change"),
      yaxis = list(title = "-log10 ajusted p-value")
    )
  
  
  volcanoPlot
    
  })
  
}




shinyApp(ui, server)