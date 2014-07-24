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
    
    checkboxInput("aggregate_on", "Aggregate antigens: ", value = FALSE),
    
    uiOutput('aggregate'),
    
    uiOutput('aggregatename'),
    
    actionButton("aggregateupdate", "Aggregate"),
    
    actionButton("updateButton","Update with filtered data"),
    
    radioButtons("method", "FDR Method: ",
                 choices = list("Expectation Maximization" = "EM", "Markov Chain Monte-Carlo" = "mcmc")
                 , selected = "EM"),
    
    downloadButton('downloadData', 'Download')
    
    
    
  ),
  
  mainPanel(    
    tabsetPanel(
      tabPanel('Data',
               dataTableOutput('data'),
               textOutput("selected"),
               plotOutput("plot"),
               uiOutput('vplotx'),
               uiOutput('vploty')
      ),
      tabPanel('FDR',
               sliderInput("threshold", "Threshold:", 
                           min = 0, max = 1, value = 0.1, step= 0.01),
               dataTableOutput('dftable'),
               radioButtons("adjustment_type", "FDR adjustment: ",
                            choices = list("Across Subjects" = "across", "Within Subjects" = "within")
                            , selected = "across"),
               plotOutput('boxplot'),
               uiOutput('xvars1'),
               uiOutput('yvars1')
               
      ),
      tabPanel('Debug',
               textOutput("debugtext")
      )
    )
    
  )
))

