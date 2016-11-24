library(ggplot2)
library(reshape2)
library(plyr)
#read in the formatted count tables
#note that I made a log table and a cpm table.  This is so that the code doesnt have to run the log calc everytime and is faster.
log2CPM <- read.table(file = "NvERTX.2_counts_plottingLOG.txt", sep="\t", header =T) #This will read in the counts table
cpm <- read.table(file = "NvERTX.2_counts_plotting.txt", sep="\t", header =T) #This will read in the counts table
Embryo <- read.table(file = "Embryo_counts_plotting_LOG.txt", sep="\t", header =T) #This will read in the counts table

#split the standard error values from the counts
log2CPMSE <- log2CPM[c(17:32)] #These split the SE data and the expression data
log2CPMCounts <- log2CPM[-c(17:32)]

colnames(log2CPMCounts) = c(-2,0,2,4,8,12,16,20,24,36,48,60,72,96,120,144) # This adds the time vector and pulls the IDS
log2CPMCounts$ID <- row.names(log2CPMCounts)
colnames(log2CPMSE) = c(-2,0,2,4,8,12,16,20,24,36,48,60,72,96,120,144)# This does the same for the standard error
log2CPMSE$ID <- row.names(log2CPMSE)

#embryo
colnames(Embryo) = c(24,48,72,96,120,144,168,192,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,7,12,24,120,240,2,7,12,24,120) # This adds the time vector and pulls the IDS
Embryo$ID <- row.names(Embryo)

#this does the whole thing over again for the cpm data
cpmSE <- cpm[c(17:32)] #These split the SE data and the expression data
cpmCounts <- cpm[-c(17:32)]

colnames(cpmCounts) = c(-2,0,2,4,8,12,16,20,24,36,48,60,72,96,120,144) # This adds the time vector and pulls the IDS
cpmCounts$ID <- row.names(cpmCounts)
colnames(cpmSE) = c(-2,0,2,4,8,12,16,20,24,36,48,60,72,96,120,144)# This does the same for the standard error
cpmSE$ID <- row.names(cpmSE)

annotations <- read.table(file = "NvERTX.2_counts_plottingAnnot.txt", sep="\t", header =T, quote = "", fill=T)
#rownames(annotations) <- annotations[,c(1)]
#annotations <- annotations[-c(1)]

server <- function(input, output, session){
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
  nve1 <- reactive({nve()[nve() != ""]})
  
  
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
    Egene()[complete.cases(Egene()),]})
  
  Egenex <- reactive({Egene1()[c(2:4, 6:10, 12:15,17:20, 22:28,34:40)]})
  
  genet <- reactive({melt(gene1())}) #this puts samples as rows, genes as columns 
  genetSE <- reactive({melt(geneSE1())}) #this puts samples as rows, genes as columns 
  Egenet <- reactive({melt(Egenex())})
  
  #this separates the Embryo datasets so they plot with different shapes
  E <- reactive({Egene1()[c(1:8,40)]})
  Et <- reactive({melt(E())})
  Fisch <- reactive({Egene1()[c(9:28,40)]})
  Ft <- reactive({melt(Fisch())})
  H <- reactive({Egene1()[c(29:34,40)]})
  Ht <- reactive({melt(H())})

  
  # this uses the SE values to make the error bar limits
  limits <- reactive(aes(ymax = genet()$value + genetSE()$value, ymin=genet()$value - genetSE()$value)) # This is the calculation for the error bars
  
  #first we output the table.  It looks nicer in the long format, hence the 't' for transpose    
  output$table <- renderTable({(gene1())})
  output$table2 <- renderTable({(annot())})
  output$table3 <- renderTable({(E())})
  output$table4 <- renderTable({(Fisch())})
  output$table5 <- renderTable({(H())})
  
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
      })
        
}
#this launches the server
#shinyApp(ui = ui, server = server)

# sweet.
