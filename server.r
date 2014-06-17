require(shiny)
require(MIMOSA)
require(ggplot2)
require(data.table)
rvalues<-reactiveValues(data=NULL,cnames=NULL,colclasses=NULL)

shinyServer(function(input, output) {
  file=reactive({
    input$file
  })

  output$data = renderDataTable({
    if(is.null(file())){
      return(NULL)
    }
    infile<-file()
    D<-fread(infile$datapath)
    rvalues$data<-D[, lapply(.SD, function(x)if(class(x)=="character")factor(x)else x),]
    
    rvalues$data
  })
  
  observe({
    rvalues$cnames<-(colnames(rvalues$data))
    })
  
  observe({
  output$colnames<-renderUI({
    cnames<-rvalues$cnames
    selectInput('colnames','Column Names',cnames)
  })
  })
  output$debugtext<-renderPrint({
    if(is.null(file())){
      return(NULL)
    }
    file()$name
  })
})
