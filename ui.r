require(shiny)
require(MIMOSA)

loadingBar <- tags$div(class="progress progress-striped active",
                       tags$div(class="bar", style="width: 100%;"))
# Code for loading message
loadingMsg <- tags$div(class="modal", tabindex="-1", role="dialog", 
                       "aria-labelledby"="myModalLabel", "aria-hidden"="true",
                       tags$div(class="modal-header",
                                tags$h3(id="myModalHeader", "Working...")),
                       tags$div(class="modal-footer",
                                loadingBar))
# The conditional panel to show when shiny is busy
loadingPanel <- conditionalPanel(paste("input.updateButton > 0 && $('html').hasClass('shiny-busy')"),
                                 loadingMsg)





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
                 , selected = "EM")
    
   
    
    
    
  ),
  
  mainPanel(    
    tabsetPanel(
      tabPanel('Data',
               dataTableOutput('data'),
               loadingPanel,
               textOutput("selected"),
               plotOutput("plot"),
               downloadLink('downloadvplot', 'Download Volcano Plot'),
               uiOutput('vplotx'),
               uiOutput('vploty')
      ),
      tabPanel('FDR',
               sliderInput("threshold", "Threshold:", 
                           min = 0, max = 0.5, value = 0.1, step= 0.01),
               dataTableOutput('dftable'),
               downloadLink('downloadfdrtable', 'Download FDR Table'),
               radioButtons("adjustment_type", "FDR adjustment: ",
                            choices = list("Across Subjects" = "across", "Within Subjects" = "within")
                            , selected = "across"),
               plotOutput('boxplot'),
               downloadLink('downloadFDRplot', 'Download FDR boxplot'),
               uiOutput('xvars1'),
               uiOutput('yvars1')
               
      ),
      tabPanel('Debug',
               textOutput("debugtext")
      )
    )
    
  )
))

