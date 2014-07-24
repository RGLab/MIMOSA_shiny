#check and load shiny, MIMOSA, ggplot2 and data.table packages
require(shiny)
require(MIMOSA)
require(ggplot2)
require(data.table)

#set rvalues to null initially
rvalues<-reactiveValues(aggregatebuttonval = 0, buttonval = 0, data=NULL,fdrtable = NULL, antigens=NULL,cytokines = NULL, visitnos = NULL, tcellsubs = NULL, xvars = NULL, yvars = NULL, threshold = 0.01, colclasses=NULL,eset=NULL,result=NULL)
rvdata_varnames = c("antigen", "cytokine", "tcellsub", "visitno")
showColumns<-c("assayid","ptid","visitno","tcellsub","cytokine","antigen","cytnum","nsub","n_antigen","cytnum_neg","nsub_neg")
#selected facet_vars
facet_vars <- list(x1 = NULL, y1 = NULL)
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
    #TAKES A SUBSET
    if(input$updateButton > rvalues$buttonval){
      if(length(input$antigens) != 0){
        rvalues$buttonval = rvalues$buttonval + 1
        thisCall <- quote(MIMOSA(data = rvalues$eset,
                                 formula = nsub+cytnum~ptid+antigen+RefTreat|cytokine+visitno+tcellsub
                                 ,method = "EM"
                                 ,subset = dummytobereassigned
                                 ,ref = dummytobereassigned))
        if(input$antigens == input$cytokines & input$visitnos == input$tcellsubs & input$antigens == input$visitnos){
          thisCall <- quote(MIMOSA(formula = nsub+cytnum~ptid+antigen+RefTreat|cytokine+visitno+tcellsub
                                   , data = rvalues$eset
                                   , method = "EM"
          ))
        }
        else if(input$antigens!= "-----" & input$cytokines != "-----" & input$visitnos != "-----"){
          thisCall[["subset"]] <- bquote(RefTreat%in%"Treatment"&antigen%in% .(input$antigens) &cytokine%in%.(input$cytokines)&visitno%in%.(input$visitnos))
          thisCall[["ref"]] <- bquote(RefTreat%in%"Reference"&antigen%in% .(input$antigens) &cytokine%in%.(input$cytokines)&visitno%in%.(input$visitnos))        
          
        }
        else if(input$antigens!= "-----" & input$cytokines == input$visitnos){
          thisCall[["subset"]] <- bquote(RefTreat%in%"Treatment"&antigen%in% .(input$antigens))
          thisCall[["ref"]] <- bquote(RefTreat%in%"Reference"&antigen%in% .(input$antigens))        
        }
        
        else if(input$cytokines != "-----" & input$antigens == input$visitnos){
          thisCall[["subset"]] <- bquote(RefTreat%in%"Treatment"&cytokine%in%.(input$cytokines))
          thisCall[["ref"]] <- bquote(RefTreat%in%"Reference"&cytokine%in%.(input$cytokines))        
        }
        
        else if(input$visitnos != "-----" & input$antigens == input$cytokines){
          thisCall[["subset"]] <- bquote(RefTreat%in%"Treatment"&visitno%in%.(input$visitnos))
          thisCall[["ref"]] <- bquote(RefTreat%in%"Reference"&visitno%in%.(input$visitnos))          
        }
        
        else if(input$antigens == "-----"){
          thisCall[["subset"]] <- bquote(RefTreat%in%"Treatment"&cytokine%in%.(input$cytokines)&visitno%in%.(input$visitnos))
          thisCall[["ref"]] <- bquote(RefTreat%in%"Reference"&cytokine%in%.(input$cytokines)&visitno%in%.(input$visitnos))            
        }
        
        else if(input$cytokines == "-----"){
          thisCall[["subset"]] <- bquote(RefTreat%in%"Treatment"&antigen%in% .(input$antigens)&visitno%in%.(input$visitnos))
          thisCall[["ref"]] <- bquote(RefTreat%in%"Reference"&antigen%in% .(input$antigens)&visitno%in%.(input$visitnos))        
        }
        
        else if(input$visitnos == "-----"){
          thisCall[["subset"]] <- bquote(RefTreat%in%"Treatment"&antigen%in% .(input$antigens) &cytokine%in%.(input$cytokines))
          thisCall[["ref"]] <- bquote(RefTreat%in%"Reference"&antigen%in% .(input$antigens) &cytokine%in%.(input$cytokines))        
        }
        
        thisCall["method"]<- bquote(.(input$method))
        #include t-cell subsetting
      }
      rvalues$result <- eval(thisCall)
    }
  })
  
  #output aggregating options if aggregating is turned on
  observe({
    if(input$aggregate_on & !is.null(rvalues$data)){
      output$aggregate<-renderUI({
        selectInput('aggregate', 'Antigen group 1', rvalues$antigens , multiple = TRUE)
      })
    }
  })
  observe({
    if(input$aggregate_on & !is.null(rvalues$data)){
      output$aggregatename<-renderUI({
        textInput('aggregatename', 'Group 1 name: ', value = input$aggregatename)
      })
    }
  })
  
  #write over chosen aggregate, delete past antigens
  observe({
    if(input$aggregate_on & !is.null(rvalues$data)){
      if(input$aggregateupdate > rvalues$aggregatebuttonval){
        rvalues$aggregatebuttonval = rvalues$aggregatebuttonval + 1
        levels(rvalues$data$antigen) <- c(levels(rvalues$data$antigen), input$aggregatename)
        rvalues$data$antigen[rvalues$data$antigen %in% input$aggregate] <- input$aggregatename
        rvalues$data <- droplevels(rvalues$data)
        rvalues$data$antigen <- factor(rvalues$data$antigen)
      }
    }
  })
  
  
  
  
  
  
  
  #generate fdrtable
  
  observe({
    if(!is.null(rvalues$result)){
      x<- rvalues$result
      q <- unlist(fdr(x), use.names = FALSE)
      pspu <- countsTable(x, proportion = TRUE)
      
      p.stim <- getZ(x)
      pd <- pData(x)
      acrossdata <-data.table(across_signif = q < (input$threshold),q, pspu, p.stim, pd)
      
      withindata<- ddply(acrossdata, .(ptid, cytokine, tcellsub, visitno), withinfunction<-function(x){
        withinframe <- data.frame(fdr_within = with(x,MIMOSA:::fdr.matrix(cbind(Pr.Nonresponse, Pr.response))))
      })
      within_signif <- withindata$fdr_within < (input$threshold)
      withindata <- data.table(within_signif, withindata)
      setkeyv(acrossdata, c("ptid", "tcellsub", "cytokine", "visitno"))
      rvalues$fdrtable <- merge(acrossdata, withindata)
      
      
      
      
    }
  })
  
  
  
  #Faceting options for VolcanoPlot
  observe({
    output$vplotx<-renderUI({
      xchoices <- setdiff(rvdata_varnames, input$vploty)
      selectInput('vplotx', 'Faceting X axis', c(".", xchoices), selected = c(".", input$vplotx), multiple = TRUE)
    })
  })
  observe({
    output$vploty<-renderUI({
      ychoices <- setdiff(rvdata_varnames, input$vplotx)
      selectInput('vploty', 'Faceting Y axis', c( ".", ychoices), selected = c(".", input$vploty), multiple = TRUE)
    })
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
  
  #display faceting options for fdr boxplot
  observe({
    output$xvars1<-renderUI({
      x1choices <- setdiff(rvdata_varnames, input$yvars1)
      selectInput('xvars1', 'Faceting X axis', c(".", x1choices), selected = c(".", input$xvars1), multiple = TRUE)
    })
  })
  observe({
    output$yvars1<-renderUI({
      y1choices <- setdiff(rvdata_varnames, input$xvars1)
      selectInput('yvars1', 'Faceting Y axis', c( ".", y1choices), selected = c(".", input$yvars1), multiple = TRUE)
    })
  })
  
  
  #printing to debug tab
  output$debugtext<-renderPrint({
    if(is.null(file())){
      return(NULL)
    }
    sprintf("%s done", input$run[1])
  })
  
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
  #rendering volcanoplot below data table
  observe({
    output$plot<-renderPlot({
      if(!is.null(rvalues$result)){
        
        if(length(input$vplotx) != 1 & length(input$vploty) == 1){
          vplotfacetformula<- as.formula(paste(". ~", addplus(input$vplotx)))
          
        }
        if(length(input$vplotx) != 1 & length(input$vploty) != 1){
          vplotfacetformula<- as.formula(paste(addplus(input$vploty), " ~", addplus(input$vplotx)))
          
          
        }
        if(length(input$vplotx)==1 & length(input$vploty) != 1){
          vplotfacetformula<- as.formula(paste(addplus(input$vploty), " ~ ."))
        }
        else if(length(input$vplotx) == 1 & length(input$vploty) == 1){
          vplotfacetformula <- NA
        }
        volcanoPlot(rvalues$result,cytnum-cytnum_REF,facet_var= vplotfacetformula)
      }
    })
    
    #display fdrtable
    output$dftable<-renderDataTable({
      if(!is.null(rvalues$result) & !is.null(rvalues$threshold)){
        rvalues$fdrtable
      }
    })
    
    #display fdr boxplot
    output$boxplot<- renderPlot({
      if(!is.null(rvalues$result)){
        plottable <- rvalues$fdrtable   
        
        if(length(input$xvars1) != 1 & length(input$yvars1) == 1){
          facetformula<- as.formula(paste(". ~", addplus(input$xvars1)))
          if(input$adjustment_type == "across"){
            p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, across_signif.fdr))  + geom_jitter(aes(colour = across_signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = across_signif.fdr)) + facet_grid(facetformula)
          }
          if(input$adjustment_type == "within"){
            p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, within_signif))  + geom_jitter(aes(colour = within_signif)) + scale_y_log10() + geom_point(data = plottable, aes(colour = within_signif)) + facet_grid(facetformula)
          }
          return(p)
        }
        if(length(input$xvars1) != 1 & length(input$yvars1) != 1){
          facetformula<- as.formula(paste(addplus(input$yvars1), " ~", addplus(input$xvars1)))
          if(input$adjustment_type == "across"){
            p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, across_signif.fdr))  + geom_jitter(aes(colour = across_signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = across_signif.fdr)) + facet_grid(facetformula)
          }
          if(input$adjustment_type == "within"){
            p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, within_signif))  + geom_jitter(aes(colour = within_signif)) + scale_y_log10() + geom_point(data = plottable, aes(colour = within_signif)) + facet_grid(facetformula)
          }
          return(p)
        }
        if(length(input$xvars1)==1 & length(input$yvars1) != 1){
          facetformula<- as.formula(paste(addplus(input$yvars1), " ~ ."))
          if(input$adjustment_type == "across"){
            p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, across_signif.fdr))  + geom_jitter(aes(colour = across_signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = across_signif.fdr)) + facet_grid(facetformula)
          }
          if(input$adjustment_type == "within"){
            p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, within_signif))  + geom_jitter(aes(colour = within_signif)) + scale_y_log10() + geom_point(data = plottable, aes(colour = within_signif)) + facet_grid(facetformula)
          }
          return(p)
        }
        
        
        if(input$adjustment_type == "across"){
          p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, across_signif.fdr))  + geom_jitter(aes(colour = across_signif.fdr)) + scale_y_log10() + geom_point(data = plottable, aes(colour = across_signif.fdr))
        }
        if(input$adjustment_type == "within"){
          p<- ggplot(plottable, aes(visitno,cytnum-cytnum_REF)) + geom_boxplot(aes(visitno, cytnum - cytnum_REF), data = subset(plottable, within_signif))  + geom_jitter(aes(colour = within_signif)) + scale_y_log10() + geom_point(data = plottable, aes(colour = within_signif))
        }
        p
      }
    })
    
    
    
    output$downloadData <- downloadHandler(
      filename = function() { paste(data, '.csv', sep='') },
      content = function(file) {
        write.csv(rvalues$fdrtable, file)
      }
    )
  })
})
