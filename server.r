library(ggplot2)
library(reshape2)
library(plyr)
library(shinyBS)
library(shinythemes)
library(Biostrings)
library(DT)
library(XML)
library(RCurl)
library(plyr)


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
  

  output$tableA = renderDataTable({
    search <- annotations[c(rows()),c(7,8,9)]
    NvERtx_ID <- row.names(search)
    datatable( cbind(' ' = '&oplus;',NvERtx_ID, search),
             rownames = T,
             escape = F,
             options = list(dom = 'ft', 
                            columnDefs = list(
                              list(visible = FALSE, targets = c(0,5)),
                              list(orderable = FALSE, className = 'details-control', targets = 1)
                            )
             ),
             callback = JS("table.column(1).nodes().to$().css({cursor: 'pointer'});
                var format = function(d) {
                  return '<div style=\"background-color:#eee; padding: .5em;\"> Extra Hits: ' +
                  d[5] + '</div>';
                  };
                table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                    row.child.hide();
                    td.html('&oplus;');
                  } else {
                    row.child(format(row.data())).show();
                    td.html('&CircleMinus;');
                  }
              });"
             ))
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    search <- annotations[c(rows()),c(7,8,9)]
    Aclicked = input$tableA_rows_selected
    updateTextInput(session, "gene1", value = row.names(search)[Aclicked[1]])
    updateTextInput(session, "gene2", value = row.names(search)[Aclicked[2]])
    updateTextInput(session, "gene3", value = row.names(search)[Aclicked[3]])
    updateTextInput(session, "gene4", value = row.names(search)[Aclicked[4]])
    updateTextInput(session, "gene5", value = row.names(search)[Aclicked[5]])
  })

  observeEvent(input$go, {
    toggleModal(session, "modal", "open")
  })
  
  #### The Mfuzz inputs
  M1 <- eventReactive(input$M1,{
    M1 <- annotations[annotations$Mfuzz_Clust %in% 1,c(4,5,7,8)] 
    M1[order(-M1$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM1 <- renderDataTable({
    M1()
  })
  
  observeEvent(input$M1, {
    toggleModal(session, "m1", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m1clicked = input$tableM1_rows_selected
      updateTextInput(session, "gene1", value = row.names(M1())[m1clicked[1]])
      updateTextInput(session, "gene2", value = row.names(M1())[m1clicked[2]])
      updateTextInput(session, "gene3", value = row.names(M1())[m1clicked[3]])
      updateTextInput(session, "gene4", value = row.names(M1())[m1clicked[4]])
      updateTextInput(session, "gene5", value = row.names(M1())[m1clicked[5]])
  })
  
  M2 <- eventReactive(input$M2,{
    M2 <- annotations[annotations$Mfuzz_Clust %in% 2,c(4,5,7,8)] 
    M2[order(-M2$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM2 <- renderDataTable({
    M2()
  })
  
  observeEvent(input$M2, {
    toggleModal(session, "m2", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m2clicked = input$tableM2_rows_selected
    updateTextInput(session, "gene1", value = row.names(M2())[m2clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M2())[m2clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M2())[m2clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M2())[m2clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M2())[m2clicked[5]])
  })
  
  M3 <- eventReactive(input$M3,{
    M3<- annotations[annotations$Mfuzz_Clust %in% 3,c(4,5,7,8)] 
    M3[order(-M3$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM3 <- renderDataTable({
    M3()
  })
  
  observeEvent(input$M3, {
    toggleModal(session, "m3", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m3clicked = input$tableM3_rows_selected
    updateTextInput(session, "gene1", value = row.names(M3())[m3clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M3())[m3clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M3())[m3clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M3())[m3clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M3())[m3clicked[5]])
  })
  
  M4 <- eventReactive(input$M4,{
    M4 <- annotations[annotations$Mfuzz_Clust %in% 4,c(4,5,7,8)] 
    M4[order(-M4$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM4 <- renderDataTable({
    M4()
  })
  
  observeEvent(input$M4, {
    toggleModal(session, "m4", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m4clicked = input$tableM4_rows_selected
    updateTextInput(session, "gene1", value = row.names(M4())[m4clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M4())[m4clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M4())[m4clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M4())[m4clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M4())[m4clicked[5]])
  })
  
  M5 <- eventReactive(input$M5,{
    M5 <- annotations[annotations$Mfuzz_Clust %in% 5,c(4,5,7,8)] 
    M5[order(-M5$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM5 <- renderDataTable({
    M5()
  })
  
  observeEvent(input$M5, {
    toggleModal(session, "m5", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m5clicked = input$tableM5_rows_selected
    updateTextInput(session, "gene1", value = row.names(M5())[m5clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M5())[m5clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M5())[m5clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M5())[m5clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M5())[m5clicked[5]])
  })
  
  M6 <- eventReactive(input$M6,{
    M6 <- annotations[annotations$Mfuzz_Clust %in% 6,c(4,5,7,8)]
    M6[order(-M6$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM6 <- renderDataTable({
    M6()
  })
  
  observeEvent(input$M6, {
    toggleModal(session, "m6", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m6clicked = input$tableM6_rows_selected
    updateTextInput(session, "gene1", value = row.names(M6())[m6clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M6())[m6clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M6())[m6clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M6())[m6clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M6())[m6clicked[5]])
  })
  
  M7 <- eventReactive(input$M7,{
    M7 <- annotations[annotations$Mfuzz_Clust %in% 7,c(4,5,7,8)]
    M7[order(-M7$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM7 <- renderDataTable({
    M7()
  })
  
  observeEvent(input$M7, {
    toggleModal(session, "m7", "open")
  })
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m7clicked = input$tableM7_rows_selected
    updateTextInput(session, "gene1", value = row.names(M7())[m7clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M7())[m7clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M7())[m7clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M7())[m7clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M7())[m7clicked[5]])
  })
  
  M8 <- eventReactive(input$M8,{
    M8 <- annotations[annotations$Mfuzz_Clust %in% 8,c(4,5,7,8)] 
    M8[order(-M8$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM8 <- renderDataTable({
    M8()
  })
  
  observeEvent(input$M8, {
    toggleModal(session, "m8", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m8clicked = input$tableM8_rows_selected
    updateTextInput(session, "gene1", value = row.names(M8())[m8clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M8())[m8clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M8())[m8clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M8())[m8clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M8())[m8clicked[5]])
  })
  
  M9 <- eventReactive(input$M9,{
    M9 <- annotations[annotations$Mfuzz_Clust %in% 9,c(4,5,7,8)]
    M9[order(-M9$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM9 <- renderDataTable({
    M9()
  })
  
  observeEvent(input$M9, {
    toggleModal(session, "m9", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m9clicked = input$tableM9_rows_selected
    updateTextInput(session, "gene1", value = row.names(M9())[m9clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M9())[m9clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M9())[m9clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M9())[m9clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M9())[m9clicked[5]])
  })
  
  M10 <- eventReactive(input$M10,{
    M10 <- annotations[annotations$Mfuzz_Clust %in% 10,c(4,5,7,8)]
    M10[order(-M10$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM10 <- renderDataTable({
    M10()
  })
  
  observeEvent(input$M10, {
    toggleModal(session, "m10", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m10clicked = input$tableM10_rows_selected
    updateTextInput(session, "gene1", value = row.names(M10())[m10clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M10())[m10clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M10())[m10clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M10())[m10clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M10())[m10clicked[5]])
  })
  
  M11 <- eventReactive(input$M11,{
    M11 <- annotations[annotations$Mfuzz_Clust %in% 11,c(4,5,7,8)]
    M11[order(-M11$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM11 <- renderDataTable({
    M11()
  })
  
  observeEvent(input$M11, {
    toggleModal(session, "m11", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m11clicked = input$tableM11_rows_selected
    updateTextInput(session, "gene1", value = row.names(M11())[m11clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M11())[m11clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M11())[m11clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M11())[m11clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M11())[m11clicked[5]])
  })
  
  M12 <- eventReactive(input$M12,{
    M12 <- annotations[annotations$Mfuzz_Clust %in% 12,c(4,5,7,8)] 
    M12[order(-M12$Mfuzz_Score),]
  },ignoreNULL= T)
  
  output$tableM12 <- renderDataTable({
    M12()
  })
  
  observeEvent(input$M12, {
    toggleModal(session, "m12", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    m12clicked = input$tableM12_rows_selected
    updateTextInput(session, "gene1", value = row.names(M12())[m12clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M12())[m12clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M12())[m12clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M12())[m12clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M12())[m12clicked[5]])
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
    c(trimws(input$gene1), trimws(input$gene2), trimws(input$gene3), trimws(input$gene4), trimws(input$gene5))
  }, ignoreNULL= T)
  
  #remove empty elements of the nve vector
  nve1 <- reactive({
    nve()[nve() != ""]
  })
  
  #subset the fasta by queries
  fasta2 <- reactive({
    a <- as.data.frame(fasta[which(names(fasta)%in%c(nve1()))])
    colnames(a) <- c('sequence')
    a
  })
  
  #output the fastas
  output$tableF <- renderTable({
    fasta2()
  },rownames = TRUE)
  

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
  #tables are reordered with Id first when printed : last column first then column from 1 to last-1
  output$table <- renderTable({
    gene1()[,c(ncol(gene1()),1:ncol(gene1())-1)]
  })
  
  # build up the annotation table with clickable links:
  annotTable <- reactive({
    uprot <- annot()$Uniprot_ID
    uprotURL <- lapply(uprot,function(i) {
      if (is.na(i)) {
        print(c("Not Available"))
      } else {
        if (i=="No_Uniprotmatch") {
          print(c("No Uniprot match"))
        } else {
        print(paste0("<a href='http://www.uniprot.org/uniprot/",i,"'>",i,"</a>"))
    }}})
    
    annotTable <- annot()
    annotTable$Uniprot_ID <- uprotURL
    
    NCBI <- strsplit(as.character(annot()$Top_NrHit..e.val), '|', fixed=TRUE)
    NCBI1 <- lapply(NCBI, '[', 1)
    NCBI2 <- lapply(NCBI, '[', 2)
    NCBI3 <- lapply(NCBI, '[', 3)
    NCBI4 <- lapply(NCBI, '[', 4)
    NCBI5 <- lapply(NCBI, '[', 5)
    NCBILink <- lapply(NCBI4,function(i) {
      if (is.na(i)) {
        print(c("No_Nr_Hit"))
      } else {
        print(paste0("<a href='https://www.ncbi.nlm.nih.gov/protein/",i,"'>",i,"</a>"))
      }})
    
    annotTable$Top_NrHit..e.val <- paste0(NCBI1,NCBI2,NCBI3,NCBILink,NCBI5)
    annotTable
  })
  
  #this makes the datatable expandle when you click the icon
  output$table2 <- DT::renderDataTable({
    datatable( cbind(' ' = '&oplus;', annotTable()),
              rownames = T,
              escape = F,
              selection = 'none',
              options = list(dom = 'ft', 
                columnDefs = list(
                  list(visible = FALSE, targets = c(0,10)),
                  list(orderable = FALSE, className = 'details-control', targets = 1)
                )
              ),
              callback = JS("table.column(1).nodes().to$().css({cursor: 'pointer'});
                var format = function(d) {
                  return '<div style=\"background-color:#eee; padding: .5em;\"> Extra Hits: ' +
                  d[10] + '</div>';
                  };
                table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                    row.child.hide();
                    td.html('&oplus;');
                  } else {
                    row.child(format(row.data())).show();
                    td.html('&CircleMinus;');
                  }
              });"
              ))
  })
  
  output$table3 <- renderTable({
    E()[,c(ncol(E()),1:ncol(E())-1)]
  })
  output$table4 <- renderTable({
    Fisch()[,c(ncol(Fisch()),1:ncol(Fisch())-1)]
  })
  output$table5 <- renderTable({
    H()[,c(ncol(H()),1:ncol(H())-1)]
  })
  
  ##output the PubMed URLs:
  pubURLs <- reactive({
    pubs <- strsplit(as.character(annot()$Top_NrHit..e.val), '|', fixed=TRUE)
    pubsp <- lapply(pubs, '[', 2)
    
    #This gets the primary reference
    urls <- lapply(pubsp,function(i) {
      if (is.na(i)) {
        print(c("Not available"))
      } else {
        print(paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed?LinkName=protein_pubmed&from_uid=",i,"'>Primary PubMed Link</a>"))
      }})
    
    #This tries to get the link to more references:
    
    more <- lapply(pubsp,function(i) {
      if (is.na(i)) {
        print(c("Not available"))
      } else {
        pubExtra <- getURL(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&dbto=pubmed&cmd=neighbor&retmode=xml&id=",i))
        pubExtraParse <- xmlParse(pubExtra)
        xml_data <- xmlToList(pubExtraParse)
        pubExtraID <- as.list(xml_data[["LinkSet"]][["LinkSetDb"]][["Link"]][["Id"]])
          if (length(pubExtraID) > 0){
            print(paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed?linkname=pubmed_pubmed&from_uid=",pubExtraID,"'>More PubMed Links</a>"))
          } else {
            print(c("Not available"))
          }
      }})
    
    urltable <- cbind(annot()$NvERTx_ID, urls, more)
    colnames(urltable) <- c("NvERTx_ID", "PubMed Link", "Extra Links")
    urltable
  })
  
  output$tableP <- DT::renderDataTable({
    datatable(pubURLs(), 
      options = list(dom = 't'))
  })
  
  ## test access col by selecting them
  output$x4 = renderPrint({
    s = input$tableA_rows_selected
  })
  ###
  
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
