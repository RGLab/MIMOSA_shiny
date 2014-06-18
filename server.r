require(shiny)
require(MIMOSA)
require(ggplot2)
require(data.table)
rvalues<-reactiveValues(data=NULL,antigens=NULL,colclasses=NULL,eset=NULL,result=NULL)
showColumns<-c("assayid","ptid","visitno","tcellsub","cytokine","antigen","cytnum","nsub","n_antigen","cytnum_neg","nsub_neg")

measure.columns=c("nsub","cytnum")
default.cast.formula=component~assayid+ptid+visitno+tcellsub+cytokine+antigen
.variable=quote(.(assayid,ptid,visitno,tcellsub,cytokine,antigen))
featureCols=1
ref.append.replace="_neg"

#MIMOSA(data = rvalues$eset,formula = nsub+cytnum~ptid+visitno+tcellsub+cytokine+antigen,method="EM",subset=RefTreat%in%"Treatment"&antigen%in%selected_antigen,ref=RefTreat%in%"Reference"&antigen%in%"selected_antigen")

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
    
    rvalues$data
  },options=list(iDisplayLength=10))
  
  observe({
    if(!is.null(rvalues$data)){
      rvalues$antigens<-levels(factor(rvalues$data$antigen))
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
  output$debugtext<-renderPrint({
    if(is.null(file())){
      return(NULL)
    }
      ant<-input$antigens
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
})
