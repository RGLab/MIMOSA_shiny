#check and load shiny, MIMOSA, ggplot2 and data.table packages
require(shiny)
require(MIMOSA)
require(ggplot2)
require(data.table)

#set rvalues to null initially
rvalues<-reactiveValues(data=NULL,antigens=NULL,cytokines = NULL, colclasses=NULL,eset=NULL,result=NULL)
#hardcoded columns?
showColumns<-c("assayid","ptid","visitno","tcellsub","cytokine","antigen","cytnum","nsub","n_antigen","cytnum_neg","nsub_neg")

measure.columns=c("nsub","cytnum")
default.cast.formula=component~assayid+ptid+visitno+tcellsub+cytokine+antigen
.variable=quote(.(assayid,ptid,visitno,tcellsub,cytokine,antigen))
featureCols=1
ref.append.replace="_neg"


shinyServer(function(input, output) {
  file=reactive({
    input$file
  })
  
  #read in file and set file to rvalues$data
  observe({
    if(is.null(file())){
      return(NULL)
    }
    infile<-file()
    D<-fread(infile$datapath)
    D<-D[, lapply(.SD, function(x)if(class(x)=="character")factor(x)else x),]
    D[,visitno:=factor(visitno)]
    D[,ptid:=factor(ptid)]
    rvalues$data<-D[,showColumns,with=FALSE]
    rvalues$eset<-ConstructMIMOSAExpressionSet(thisdata=data.frame(D),
                                               reference = NULL,
                                               measure.columns = measure.columns,
                                               default.cast.formula=default.cast.formula,
                                               .variables = .variable,
                                               featureCols=featureCols,
                                               ref.append.replace = ref.append.replace)
  })
  
  output$data = renderDataTable({
    if(is.null(rvalues$data)){
      return(NULL)
    }
    
    #filter rvalues$data for antigens and cytokines
    if(is.null(input$antigens)&is.null(input$cytokines)){
      return(rvalues$data)
    }
    if(is.null(input$antigens)){
      return(rvalues$data[rvalues$data$cytokine == input$cytokines])
    }
    if(is.null(input$cytokines)){
      return(rvalues$data[rvalues$data$antigen == input$antigens])
    }
    rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$cytokine == input$cytokines ]
    
    
  }
  ,options=list(iDisplayLength=10)
  )
  
  observe({
    if(!is.null(rvalues$data)){
      rvalues$antigens<-levels(factor(rvalues$data$antigen))
      rvalues$cytokines<-levels(factor(rvalues$data$cytokine))
    }
  })
  
  observe({
    output$antigens<-renderUI({
      if(is.null(rvalues$antigens)){
        return(NULL)
      }
      antigens<-rvalues$antigens
      selectInput('antigens','Antigens',antigens)
    })
  })
  observe({
    output$cytokines<-renderUI({
      if(is.null(rvalues$cytokines)){
        return(NULL)
      }
      cytokines<-rvalues$cytokines
      selectInput('cytokines','Cytokines',cytokines)
    })
  })
  
  
  output$debugtext<-renderPrint({
    if(is.null(file())){
      return(NULL)
    }
    rvalues$result<-MIMOSA(data = rvalues$eset,formula = nsub+cytnum~ptid+antigen+RefTreat|cytokine+visitno+tcellsub,method="EM",subset=RefTreat%in%"Treatment"&antigen%in%"CMV"&cytokine%in%"IL2+"&visitno%in%"10",ref=RefTreat%in%"Reference"&antigen%in%"CMV"&cytokine%in%"IL2+"&visitno%in%"10")
    sprintf("%s done", input$run[1])
  })
  
  observe({
    output$plot<-renderPlot({
      if(!is.null(rvalues$result)){
        volcanoPlot(rvalues$result,cytnum-cytnum_REF,facet_var=~tcellsub)
      }
    })
  })
