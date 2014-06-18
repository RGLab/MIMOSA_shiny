require(shiny)
require(MIMOSA)
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("MIMOSA: Mixture Models for Single Cell Assays"),
  
  sidebarPanel(
    fileInput('file','Choose Text File in Standard SCHARP Format',accept = c('text/csv','.csv')),
    uiOutput('antigens'),
    actionButton("run","Run MIMOSA")
    
    ),
  
  mainPanel(    
    tabsetPanel(
        tabPanel('Data',
          dataTableOutput('data'),
          plotOutput("plot")
          ),
        tabPanel('Debug',
          textOutput("debugtext")
        )
      )
      
  )
))
