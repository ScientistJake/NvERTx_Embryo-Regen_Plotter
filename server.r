library(ggplot2)
library(reshape2)
library(plyr)
library(shinyBS)

server <- function(input, output, session){
  
  

  withProgress(message = 'Loading Data',
               detail = 'This may take a while...', value = 0, {
                 load(file = "Data.RData")
              })
  
  rows <- eventReactive(input$go,{
    which(
      rowSums(
        `dim<-`(grepl(input$search, as.matrix(annotations), ignore.case=TRUE), dim(annotations))
      ) > 0
    )
  },ignoreNULL= T)

  
  output$tableA <- renderTable({
    annotations[c(rows()),c(1,7,8,9)]
  })
  
  
  observeEvent(input$go, {
    toggleModal(session, "modal", "open")
  })
  
  #### The Mfuzz inputs
   M1 <- eventReactive(input$M1,{
      annotations[annotations$Mfuzz_Clust %in% 1,c(1,4,5,7,8)] 
  },ignoreNULL= T)
   
   output$tableM1 <- renderTable({
     M1()[order(-M1()$Mfuzz_Score),]
   })
   
  observeEvent(input$M1, {
    toggleModal(session, "m1", "open")
  })
  
  M2 <- eventReactive(input$M2,{
    annotations[annotations$Mfuzz_Clust %in% 2,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM2 <- renderTable({
    M2()[order(-M2()$Mfuzz_Score),]
  })
  
  observeEvent(input$M2, {
    toggleModal(session, "m2", "open")
  })
  
  M3 <- eventReactive(input$M3,{
    annotations[annotations$Mfuzz_Clust %in% 3,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM3 <- renderTable({
    M3()[order(-M3()$Mfuzz_Score),]
  })
  
  observeEvent(input$M3, {
    toggleModal(session, "m3", "open")
  })
  
  M4 <- eventReactive(input$M4,{
    annotations[annotations$Mfuzz_Clust %in% 4,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM4 <- renderTable({
    M4()[order(-M4()$Mfuzz_Score),]
  })
  
  observeEvent(input$M4, {
    toggleModal(session, "m4", "open")
  })
  
  M5 <- eventReactive(input$M5,{
    annotations[annotations$Mfuzz_Clust %in% 5,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM5 <- renderTable({
    M5()[order(-M5()$Mfuzz_Score),]
  })
  
  observeEvent(input$M5, {
    toggleModal(session, "m5", "open")
  })
  
  M6 <- eventReactive(input$M6,{
    annotations[annotations$Mfuzz_Clust %in% 6,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM6 <- renderTable({
    M6()[order(-M6()$Mfuzz_Score),]
  })
  
  observeEvent(input$M6, {
    toggleModal(session, "m6", "open")
  })
  
  M7 <- eventReactive(input$M7,{
    annotations[annotations$Mfuzz_Clust %in% 7,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM7 <- renderTable({
    M7()[order(-M7()$Mfuzz_Score),]
  })
  
  observeEvent(input$M7, {
    toggleModal(session, "m7", "open")
  })
  
  M8 <- eventReactive(input$M8,{
    annotations[annotations$Mfuzz_Clust %in% 8,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM8 <- renderTable({
    M8()[order(-M8()$Mfuzz_Score),]
  })
  
  observeEvent(input$M8, {
    toggleModal(session, "m8", "open")
  })
  
  M9 <- eventReactive(input$M9,{
    annotations[annotations$Mfuzz_Clust %in% 9,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM9 <- renderTable({
    M9()[order(-M9()$Mfuzz_Score),]
  })
  
  observeEvent(input$M9, {
    toggleModal(session, "m9", "open")
  })
  
  M10 <- eventReactive(input$M10,{
    annotations[annotations$Mfuzz_Clust %in% 10,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM10 <- renderTable({
    M10()[order(-M10()$Mfuzz_Score),]
  })
  
  observeEvent(input$M10, {
    toggleModal(session, "m10", "open")
  })
  
  M11 <- eventReactive(input$M11,{
    annotations[annotations$Mfuzz_Clust %in% 11,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM11 <- renderTable({
    M11()[order(-M11()$Mfuzz_Score),]
  })
  
  observeEvent(input$M11, {
    toggleModal(session, "m11", "open")
  })
  
  M12 <- eventReactive(input$M12,{
    annotations[annotations$Mfuzz_Clust %in% 12,c(1,4,5,7,8)] 
  },ignoreNULL= T)
  
  output$tableM12 <- renderTable({
    M12()[order(-M12()$Mfuzz_Score),]
  })
  
  observeEvent(input$M12, {
    toggleModal(session, "m12", "open")
  })
  
  #End of Mfuzz inputs
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['gene1']])) {
      updateTextInput(session, "gene1", value = query[['gene1']])
    }
  })
  #this responds to the checkbox and decides whether to use the log data or cpm
  counts <- reactive({
    if (as.numeric(input$log) >0 ){
      return(log2CPMCounts)
    } else {
      return(cpmCounts)
    }
  })
  
  #same for the standard error values
  countsSE <- reactive({
    if (as.numeric(input$log) >0 ){
      return(log2CPMSE)
    } else {
      return(cpmSE)
    }
  })
  
  #and we want to change the axis label for the plot if log is ticked
  y_label <- reactive({
    if (as.numeric(input$log) >0 ){
      return(c("Log2(Counts per million +1)"))
    } else {
      return(c("Counts per million"))
    }
  })
  

  #the event reactive call here makes the code wait until the user clicks the evaluate button.
  #this assigns the NvERTx numbers to the nve object
  nve <- eventReactive(input$do, {
    c(input$gene1, input$gene2, input$gene3, input$gene4, input$gene5)
  }, ignoreNULL= T)
  
  #remove empty elements of the nve vector
  nve1 <- reactive({
    nve()[nve() != ""]
  })
  
  
  #this subsets the count table by the nve numbers
  gene <- reactive({
    counts()[c(nve1()),]
  })
  
  #embryo
  Egene <- reactive({
    Embryo[c(nve1()),]
  })
  
  
  #subsetting the annotations
  annot <- reactive({
    annotations[c(nve1()),]
  })
  #this removes all NA, in case the users puts in a number with no hits.
  
  
  #note that I am NOT overwriting gene. I'm passing it to gene1. This avoids some crazy infinite loop error.
  
  
  gene1 <- reactive({
    gene()[complete.cases(gene()),]
  })
  
  #same workflow, subsetting the SE values
  geneSE <- reactive({
    countsSE()[c(nve1()),]
  })
  
  #removing NAs
  geneSE1 <- reactive({
    geneSE()[complete.cases(geneSE()),]
  })
  
  #Extra step for embryo to get rid of duplicate values
  Egene1 <- reactive({
    Egene()[complete.cases(Egene()),]
  })
  
  Egenex <- reactive({
    Egene1()[c(2:4, 6:10, 12:15,17:20, 22:28,34:40)]
  })
  
  genet <- reactive({
    melt(gene1())
  }) #this puts samples as rows, genes as columns 
  genetSE <- reactive({
    melt(geneSE1())
  }) #this puts samples as rows, genes as columns 
  Egenet <- reactive({
    melt(Egenex())
  })
  
  #this separates the Embryo datasets so they plot with different shapes
  E <- reactive({
    Egene1()[c(1:8,40)]
  })
  Et <- reactive({
    melt(E())
  })
  Fisch <- reactive({
    Egene1()[c(9:28,40)]
  })
  Ft <- reactive({
    melt(Fisch())
  })
  H <- reactive({
    Egene1()[c(29:34,40)]
  })
  Ht <- reactive({
    melt(H())
  })

  
  # this uses the SE values to make the error bar limits
  limits <- reactive(aes(ymax = genet()$value + genetSE()$value, ymin=genet()$value - genetSE()$value)) # This is the calculation for the error bars
  
  #first we output the table.  It looks nicer in the long format, hence the 't' for transpose    
  output$table <- renderTable({
    (gene1())
  })
  output$table2 <- renderTable({
    annot()
  })
  output$table3 <- renderTable({
    E()
  })
  output$table4 <- renderTable({
    Fisch()
  })
  output$table5 <- renderTable({
    H()
  })
  
  #now we output the plot.  
  #this is the regen:
  p <- reactive({
          ggplot(genet(), aes(x=as.numeric(as.character(genet()$variable)), y=value, colour=ID)) + geom_line() +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
     
  })
 
  output$plot1 <- renderPlot({ 
    #to get the graph to show up in shiny you need to print 
    print(p()+ scale_x_continuous(minor_breaks = NULL, breaks=c(0,2,4,8,12,16,20,24,36,48,60,72,96,120,144)) +
            geom_errorbar(limits(), width=0.2) +
            ylab(y_label()) +
            xlab("Hours Post Amputation"))  
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'Rplot.pdf',
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::pdf(..., width = as.numeric(input$w), height=as.numeric(input$h)
                       )
      }
      ggsave(file, plot = p()+ scale_x_continuous(minor_breaks = NULL, breaks=c(0,2,4,8,12,16,20,24,36,48,60,72,96,120,144)) +
               geom_errorbar(limits(), width=0.2) +
               ylab(y_label()) +
               xlab("Hours Post Amputation"), device = device)
  })
  #embryo plot
  q <- reactive({
    ggplot(Egenet(), aes(x=as.numeric(as.character(Egenet()$variable)), y=value, colour=ID)) + geom_line() +
      scale_shape_discrete(solid=F, name='Dataset') +
      geom_point(data = Et(), aes(x=as.numeric(as.character(Et()$variable)), y=value, colour=ID, shape='Warner et al. (2017)'), size=2) +
      geom_point(data = Ft(), aes(x=as.numeric(as.character(Ft()$variable)), y=value, colour=ID, shape='Fischer et al. (2014)'),size=2) +
      geom_point(data = Ht(), aes(x=as.numeric(as.character(Ht()$variable)), y=value, colour=ID, shape='Helm et al. (2103)'),size=2) +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
  })
  output$plot2 <- renderPlot({ 
      #to get the graph to show up in shiny you need to print 
      print(q()+ scale_x_continuous(minor_breaks = NULL, breaks=c(0,6,12,24,48,72,96,120,144,168,192,120,240)) +
              ylab("Log2 CPM") +
              xlab("Hours Post Fertilization"))  
  })
  output$downloadPlot2 <- downloadHandler(
      filename = 'Eplot.pdf',
      content = function(file) {
        device <- function(..., width, height) {
          grDevices::pdf(..., width = as.numeric(input$w), height=as.numeric(input$h)
          )
        }
        ggsave(file, plot = q()+ scale_x_continuous(minor_breaks = NULL, breaks=c(0,6,12,24,48,72,96,120,144,168,192,120,240)) +
                 ylab("Log2(Counts per million +1)") +
                 xlab("Hours Post Fertilization"), device = device)
      }
  )
        
}


#this launches the server
#shinyApp(ui = ui, server = server)

# sweet.
