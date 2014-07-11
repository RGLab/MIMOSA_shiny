#check and load shiny, MIMOSA, ggplot2 and data.table packages
require(shiny)
require(MIMOSA)
require(ggplot2)
require(data.table)

#set rvalues to null initially
rvalues<-reactiveValues(data=NULL,antigens=NULL,cytokines = NULL, visitnos = NULL, tcellsubs = NULL, xvars = NULL, yvars = NULL, threshold = 0.01, colclasses=NULL,eset=NULL,result=NULL)
rvdata_varnames = c("antigen", "cytokine", "tcellsub", "visitno")
#hardcoded columns?
showColumns<-c("assayid","ptid","visitno","tcellsub","cytokine","antigen","cytnum","nsub","n_antigen","cytnum_neg","nsub_neg")
#selected facet_vars
facet_vars <- list(x1 = NULL, y1 = NULL, x2 = NULL, y2 = NULL)
measure.columns=c("nsub","cytnum")
default.cast.formula=component~assayid+ptid+visitno+tcellsub+cytokine+antigen
.variable=quote(.(assayid,ptid,visitno,tcellsub,cytokine,antigen))
featureCols=1
ref.append.replace="_neg"


shinyServer(function(input, output) {
  file=reactive({
    input$file
  })
  
  #read in file and set file to rvalues$data, pass to MIMOSA expression set.
  #store result values to be used in plot
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
    thisdata = data.frame(D)
    rvalues$eset<-ConstructMIMOSAExpressionSet(thisdata,
                                               reference = NULL,
                                               measure.columns = measure.columns,
                                               default.cast.formula=default.cast.formula,
                                               .variables = .variable,
                                               featureCols=featureCols,
                                               ref.append.replace = ref.append.replace)
    rvalues$result<-MIMOSA(data = rvalues$eset,formula = nsub+cytnum~ptid+antigen+RefTreat|cytokine+visitno+tcellsub,method="EM",subset=RefTreat%in%"Treatment"&antigen%in%"CMV"&cytokine%in%"IL2+"&visitno%in%"10",ref=RefTreat%in%"Reference"&antigen%in%"CMV"&cytokine%in%"IL2+"&visitno%in%"10")
  })
  
  #filter data based on dropdown input and print to table
  output$data = renderDataTable({
    if(is.null(rvalues$data)){
      return(NULL)
    }
    
    #filter rvalues$data for dropdown input
    #terrible ugly looking chunk of if statements
    if(is.null(input$antigens)|is.null(input$cytokines)|is.null(input$visitnos)|is.null(input$tcellsubs)){
      return(rvalues$data)
    }
    
    #if all dropdowns are set to none
    if(input$antigens == "-----" & input$cytokines == "-----" & input$visitnos == "-----" & input$tcellsubs =="-----"){
      return(rvalues$data)
    }
    
    #if three dropdowns are set to none
    else if(input$antigens == "-----" & input$cytokines == "-----" & input$visitnos == "-----"){
      return(rvalues$data[rvalues$data$tcellsub == input$tcellsubs])
    }
    else if(input$antigens == "-----" & input$cytokines == "-----"& input$tcellsubs =="-----"){
      return(rvalues$data[rvalues$data$visitno == input$visitnos])
    }
    else if(input$antigens == "-----" & input$visitnos == "-----" & input$tcellsubs =="-----"){
      return(rvalues$data[rvalues$data$cytokine == input$cytokines])
    }
    else if(input$cytokines == "-----" & input$visitnos == "-----" & input$tcellsubs =="-----"){
      return(rvalues$data[rvalues$data$antigen == input$antigens])
    }
    
    #if two dropdowns are set to none
    else if(input$antigens == "-----" & input$cytokines =="-----"){
      return(rvalues$data[rvalues$data$visitno == input$visitnos & rvalues$data$tcellsub == input$tcellsubs])
    }
    else if(input$antigens == "-----" & input$visitnos =="-----"){
      return(rvalues$data[rvalues$data$cytokine == input$cytokines & rvalues$data$tcellsub == input$tcellsubs])
    }
    else if(input$antigens == "-----" & input$tcellsubs =="-----"){
      return(rvalues$data[rvalues$data$cytokine == input$cytokines & rvalues$data$visitno == input$visitnos])
    }
    else if(input$cytokines == "-----" & input$visitnos =="-----"){
      return(rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$tcellsub == input$tcellsubs])
    }
    else if(input$cytokines == "-----" & input$tcellsubs =="-----"){
      return(rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$visitno == input$visitnos])
    }
    else if(input$visitnos == "-----" & input$tcellsubs =="-----"){
      return(rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$cytokine == input$cytokines])
    }
    
    #if one set to none
    else if(input$antigens == "-----"){
      return(rvalues$data[rvalues$data$cytokine == input$cytokines & rvalues$data$visitno == input$visitnos & rvalues$data$tcellsub == input$tcellsubs])
    }
    else if(input$cytokines == "-----"){
      return(rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$visitno == input$visitnos & rvalues$data$tcellsub == input$tcellsubs])
    }
    else if(input$visitnos == "-----"){
      return(rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$cytokine == input$cytokines
                          & rvalues$data$tcellsub == input$tcellsubs])
    }
    else if(input$tcellsubs == "-----"){
      return(rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$cytokine == input$cytokines
                          & rvalues$data$visitno == input$visitnos])
    }
    rvalues$data[rvalues$data$antigen == input$antigens & rvalues$data$cytokine == input$cytokines & rvalues$data$visitno == input$visitnos & rvalues$data$tcellsub == input$tcellsubs]
  }
  ,options=list(iDisplayLength=10)
  )
  
  #? clean up values: antigen, cytokine, visitno, tcellsub
  observe({
    if(!is.null(rvalues$data)){
      rvalues$antigens<-levels(factor(rvalues$data$antigen))
      rvalues$cytokines<-levels(factor(rvalues$data$cytokine))
      rvalues$visitnos<-levels(factor(rvalues$data$visitno))
      rvalues$tcellsubs<-levels(factor(rvalues$data$tcellsub))
      
    }
  })
  
  #set dropdowns to print when data is in
  observe({
    output$antigens<-renderUI({
      if(is.null(rvalues$antigens)){
        return(NULL)
      }
      antigens<-rvalues$antigens
      selectInput('antigens','Antigens',c("-----", antigens))
    })
    output$cytokines<-renderUI({
      if(is.null(rvalues$cytokines)){
        return(NULL)
      }
      cytokines<-rvalues$cytokines
      selectInput('cytokines','Cytokines',c("-----",cytokines))
    })
    output$visitnos<-renderUI({
      if(is.null(rvalues$visitnos)){
        return(NULL)
      }
      visitnos<-rvalues$visitnos
      selectInput('visitnos','Visit Numbers',c("-----",visitnos))
    })
    output$tcellsubs<-renderUI({
      if(is.null(rvalues$tcellsubs)){
        return(NULL)
      }
      tcellsubs<-rvalues$tcellsubs
      selectInput('tcellsubs','T cell subs',c("-----",tcellsubs))
    })
  })
  
  
  observe({
    output$xvars1<-renderUI({
      x1choices <- setdiff(rvdata_varnames, input$yvars1)
      selectInput('xvars1', 'X variable', c(".", x1choices), selected = c(".", input$xvars1), multiple = TRUE)
    })
  })
  observe({
    output$yvars1<-renderUI({
      y1choices <- setdiff(rvdata_varnames, input$xvars1)
      selectInput('yvars1', 'Y variable', c( ".", y1choices), selected = c(".", input$yvars1), multiple = TRUE)
    })
  })
  
  
  #printing to debug tab
  output$debugtext<-renderPrint({
    if(is.null(file())){
      return(NULL)
    }
    sprintf("%s done", input$run[1])
  })
  
  #rendering volcanoplot below data table
  observe({
    output$plot<-renderPlot({
      if(!is.null(rvalues$result)){
        volcanoPlot(rvalues$result,cytnum-cytnum_REF,facet_var=~tcellsub)
      }
    })
    
    output$dftable<-renderDataTable({
      if(!is.null(rvalues$result) & !is.null(rvalues$threshold)){
        x<- rvalues$result
        q <- -log10(unlist(fdr(x), use.names = FALSE))
        pspu <- countsTable(x, proportion = TRUE)
        p.stim <- getZ(x)
        pd <- pData(x)
        fdrtable <-data.table(signif = q > -log10(input$threshold),q, pspu, p.stim, pd)
        fdrtable
      }
    })
    
    output$boxplot<- renderPlot({
      if(!is.null(rvalues$result)){
        x<- rvalues$result
        q <- -log10(unlist(fdr(x), use.names = FALSE))
        pspu <- countsTable(x, proportion = TRUE)
        p.stim <- getZ(x)
        pd <- pData(x)
        plottable <-data.table(signif = q > -log10(input$threshold),q, pspu, p.stim, pd)
        
        
        
        #ADDPLUS function here
        addplus <- function(input){
          size<- length(input)
          if(size== 1){
            return(input)
          }
          retval<- input[1]
          for(i in 2:size){
            retval <- paste(retval, " + ", input[i])
          }
          return(retval)
        }
        
        
          if(length(input$xvars1) != 1 & length(input$yvars1) == 1){
          facetformula<- as.formula(paste(". ~", addplus(input$xvars1)))
          
          p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, signif.fdr))  + geom_jitter(aes(colour = signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = signif.fdr)) + facet_grid(facetformula)
          return(p)
          }
        if(length(input$xvars1) != 1 & length(input$yvars1) != 1){
          facetformula<- as.formula(paste(addplus(input$yvars1), " ~", addplus(input$xvars1)))
          
          p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, signif.fdr))  + geom_jitter(aes(colour = signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = signif.fdr)) + facet_grid(facetformula)
          return(p)
        }
        if(length(input$xvars1)==1 & length(input$yvars1) != 1){
          facetformula<- as.formula(paste(addplus(input$yvars1), " ~ ."))
          
          p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, signif.fdr))  + geom_jitter(aes(colour = signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = signif.fdr)) + facet_grid(facetformula)
          return(p)
        }
          
          
        p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, signif.fdr))  + geom_jitter(aes(colour = signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = signif.fdr))
        p
      }
    })
  })
})
