require(shiny)
require(MIMOSA)
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("MIMOSA: Mixture Models for Single Cell Assays"),
  
  sidebarPanel(
    fileInput('file','Choose Text File in Standard SCHARP Format',accept = c('text/csv','.csv')),
    uiOutput('colnames')
    ),
  
  mainPanel(    
    tabsetPanel(
        tabPanel('Data',
          dataTableOutput('data')
          ),
        tabPanel('Debug',
          textOutput("debugtext")
        )
      )
      
  )
))
