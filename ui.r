require(shiny)
require(MIMOSA)
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("MIMOSA: Mixture Models for Single Cell Assays"),
  
  sidebarPanel(
    fileInput('file','Choose Text File in Standard SCHARP Format',accept = c('text/csv','.csv')),
    uiOutput('antigens'),
    uiOutput('cytokines'),
    uiOutput('visitnos'),
    uiOutput('tcellsubs'),
    
    actionButton("run","Run MIMOSA")
    
  ),
  
  mainPanel(    
    tabsetPanel(
      tabPanel('Data',
               dataTableOutput('data'),
               #
               textOutput("selected"),
               plotOutput("plot")
      ),
      tabPanel('Debug',
               textOutput("debugtext")
      )
    )
    
  )
))

