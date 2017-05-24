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
library(emoGG)


server <- function(input, output, session){
  
#  output$frame <- renderUI({
#    my_test <- tags$iframe(src="http://134.59.51.195:4567/", height=600, width=1000)
#    print(my_test)
#    my_test
#  })
  
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
  

  #This is the converstion Tool for NvERTX.2 to NvERTx.4

  NewIDS <- eventReactive(input$convert,{
    old_IDS <- trimws(unlist(strsplit(as.character(input$NveConvert), ',', fixed=TRUE)))
    IDtable[IDtable$NvERTx.2 %in% old_IDS,]
  },ignoreNULL= T)
  
  output$conversionTable <- renderDataTable({
    NewIDS()
  }, rownames= FALSE)
  
  observeEvent(input$convert, {
    toggleModal(session, "converted", "open")
  })
  
  observe({
    NewIDSclicked = input$conversionTable_rows_selected
    updateTextInput(session, "gene1", value = NewIDS()$NvERTx.4[NewIDSclicked[1]])
    updateTextInput(session, "gene2", value = NewIDS()$NvERTx.4[NewIDSclicked[2]])
    updateTextInput(session, "gene3", value = NewIDS()$NvERTx.4[NewIDSclicked[3]])
    updateTextInput(session, "gene4", value = NewIDS()$NvERTx.4[NewIDSclicked[4]])
    updateTextInput(session, "gene5", value = NewIDS()$NvERTx.4[NewIDSclicked[5]])
  })
  
  output$tableA = renderDataTable({
    search <- annotations[c(rows()),c(9,10,11)]
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
    M1 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-1"),c(4,5,9,10)] 
    M1[order(-M1$Mfuzz_R_Score),]
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
    M2 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-2"),c(4,5,9,10)] 
    M2[order(-M2$Mfuzz_R_Score),]
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
    M3<- annotations[annotations$Mfuzz_R_Clust %in% c("R-3"),c(4,5,9,10)] 
    M3[order(-M3$Mfuzz_R_Score),]
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
    M4 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-4"),c(4,5,9,10)] 
    M4[order(-M4$Mfuzz_R_Score),]
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
    M5 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-5"),c(4,5,9,10)] 
    M5[order(-M5$Mfuzz_R_Score),]
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
    M6 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-6"),c(4,5,9,10)]
    M6[order(-M6$Mfuzz_R_Score),]
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
    M7 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-7"),c(4,5,9,10)]
    M7[order(-M7$Mfuzz_R_Score),]
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
    M8 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-8"),c(4,5,9,10)] 
    M8[order(-M8$Mfuzz_R_Score),]
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
    M9 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-9"),c(4,5,9,10)]
    M9[order(-M9$Mfuzz_R_Score),]
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
    M10 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-10"),c(4,5,9,10)]
    M10[order(-M10$Mfuzz_R_Score),]
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

  #
  #
  # Embryonic Mfuzz
  #
  #
  
  E1 <- eventReactive(input$E1,{
    E1 <- annotations[annotations$Mfuzz_E_Clust %in% c("E-1"),c(6,7,9,10)] 
    E1[order(-E1$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE1 <- renderDataTable({
    E1()
  })
  
  observeEvent(input$E1, {
    toggleModal(session, "e1", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e1clicked = input$tableE1_rows_selected
    updateTextInput(session, "gene1", value = row.names(E1())[e1clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E1())[e1clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E1())[e1clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E1())[e1clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E1())[e1clicked[5]])
  })
  
  E2 <- eventReactive(input$E2,{
    E2 <- annotations[annotations$Mfuzz_E_Clust %in% c("E-2"),c(6,7,9,10)] 
    E2[order(-E2$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE2 <- renderDataTable({
    E2()
  })
  
  observeEvent(input$E2, {
    toggleModal(session, "e2", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e2clicked = input$tableE2_rows_selected
    updateTextInput(session, "gene1", value = row.names(E2())[e2clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E2())[e2clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E2())[e2clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E2())[e2clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E2())[e2clicked[5]])
  })
  
  E3 <- eventReactive(input$E3,{
    E3<- annotations[annotations$Mfuzz_E_Clust %in% c("E-3"),c(6,7,9,10)] 
    E3[order(-E3$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE3 <- renderDataTable({
    E3()
  })
  
  observeEvent(input$E3, {
    toggleModal(session, "e3", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e3clicked = input$tableM3_rows_selected
    updateTextInput(session, "gene1", value = row.names(E3())[e3clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E3())[e3clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E3())[e3clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E3())[e3clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E3())[e3clicked[5]])
  })
  
  E4 <- eventReactive(input$E4,{
    E4 <- annotations[annotations$Mfuzz_E_Clust %in% c("E-4"),c(6,7,9,10)] 
    E4[order(-E4$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE4 <- renderDataTable({
    E4()
  })
  
  observeEvent(input$E4, {
    toggleModal(session, "e4", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e4clicked = input$tableE4_rows_selected
    updateTextInput(session, "gene1", value = row.names(E4())[e4clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E4())[e4clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E4())[e4clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E4())[e4clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E4())[e4clicked[5]])
  })
  
  E5 <- eventReactive(input$E5,{
    E5 <- annotations[annotations$Mfuzz_E_Clust %in% c("E-5"),c(6,7,9,10)] 
    E5[order(-E5$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE5 <- renderDataTable({
    E5()
  })
  
  observeEvent(input$E5, {
    toggleModal(session, "e5", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e5clicked = input$tableE5_rows_selected
    updateTextInput(session, "gene1", value = row.names(E5())[e5clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E5())[e5clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E5())[e5clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E5())[e5clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E5())[e5clicked[5]])
  })
  
  E6 <- eventReactive(input$E6,{
    E6 <- annotations[annotations$Mfuzz_E_Clust %in% c("E-6"),c(6,7,9,10)]
    E6[order(-E6$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE6 <- renderDataTable({
    E6()
  })
  
  observeEvent(input$E6, {
    toggleModal(session, "e6", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e6clicked = input$tableE6_rows_selected
    updateTextInput(session, "gene1", value = row.names(E6())[e6clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E6())[e6clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E6())[e6clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E6())[e6clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E6())[e6clicked[5]])
  })
  
  E7 <- eventReactive(input$E7,{
    E7 <- annotations[annotations$Mfuzz_E_Clust %in% c("E-7"),c(6,7,9,10)]
    E7[order(-E7$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE7 <- renderDataTable({
    E7()
  })
  
  observeEvent(input$E7, {
    toggleModal(session, "e7", "open")
  })
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e7clicked = input$tableE7_rows_selected
    updateTextInput(session, "gene1", value = row.names(E7())[e7clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E7())[e7clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E7())[e7clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E7())[e7clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E7())[e7clicked[5]])
  })
  
  E8 <- eventReactive(input$E8,{
    E8 <- annotations[annotations$Mfuzz_E_Clust %in% c("E-8"),c(6,7,9,10)] 
    E8[order(-E8$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE8 <- renderDataTable({
    E8()
  })
  
  observeEvent(input$E8, {
    toggleModal(session, "e8", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    e8clicked = input$tableE8_rows_selected
    updateTextInput(session, "gene1", value = row.names(E8())[e8clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E8())[e8clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E8())[e8clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E8())[e8clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E8())[e8clicked[5]])
  })
  #
  #
  #
  #
  #
  #End of Mfuzz inputs
  #
  #
  #
  
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
  #same workflow, subsetting the SE values
  EgeneSE <- reactive({
    EmbryoSE[c(nve1()),]
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
  ERSE <- reactive({
    ERSEx <- EgeneSE()[c(1:8,35)]
    melt(ERSEx)
  })
  
  Fisch <- reactive({
    Egene1()[c(9:28,40)]
  })
  Ft <- reactive({
    melt(Fisch())
  })
  FSE <- reactive({
    FSEx <- EgeneSE()[c(9:28,35)]
    melt(FSEx)
  })
  
  
  H <- reactive({
    Egene1()[c(29:34,40)]
  })
  Ht <- reactive({
    melt(H())
  })
  HSE <- reactive({
    HSEx <- EgeneSE()[c(29:34,35)]
    melt(HSEx)
  })
  
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
                  list(visible = FALSE, targets = c(0,12)),
                  list(orderable = FALSE, className = 'details-control', targets = 1)
                )
              ),
              callback = JS("table.column(1).nodes().to$().css({cursor: 'pointer'});
                var format = function(d) {
                  return '<div style=\"background-color:#eee; padding: .5em;\"> Extra Hits: ' +
                  d[12] + '</div>';
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
  
  
  ##
  ## Plots
  ##
  ##
  
  # this uses the SE values to make the error bar limits
  limits <- reactive(aes(ymax = genet()$value + genetSE()$value, ymin=genet()$value - genetSE()$value)) # This is the calculation for the error bars
  #for whatever reason I couldnt do the same for the embryonic. I had to explicitly map them in the Geom_errorbar() below
  
  
  #now we output the plot.  
  #this is the regen:
  p <- reactive({
    ggplot(genet(), aes(x=as.numeric(as.character(genet()$variable)), y=value, colour=ID)) + geom_line() +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
    
  })
  r <- reactive({
    ggplot(genet(), aes(x=as.numeric(as.character(genet()$variable)), y=value, colour=ID)) + geom_line() +
      geom_emoji(emoji="1f363") +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
  })
  s <- reactive({
    ggplot(genet(), aes(x=as.numeric(as.character(genet()$variable)), y=value, colour=ID)) + geom_line() +
      geom_emoji(emoji="1f431") +
    theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
  })
  t <- reactive({
    ggplot(genet(), aes(x=as.numeric(as.character(genet()$variable)), y=value, colour=ID)) + geom_line() +
      geom_emoji(emoji="1f419") +
    theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
          axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
  })
  
  rplot <- reactive({
    emog <- sample(c(1:100), 1, replace = FALSE, prob = NULL)
    if (emog %in% c(1)){
      return(r())
    } else if (emog %in% c(2)){
      return(s())
    } else if (emog %in% c(3)){
      return(t())
    } else {
      return(p())
    }
  })
  
  
  output$plot1 <- renderPlot({ 
    #to get the graph to show up in shiny you need to print 
    print(rplot()+ scale_x_continuous(minor_breaks = NULL, breaks=c(0,2,4,8,12,16,20,24,36,48,60,72,96,120,144)) +
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
            geom_errorbar(data = Et(), aes(x=as.numeric(as.character(Et()$variable)), y=value, ymax = Et()$value + ERSE()$value, ymin=Et()$value - ERSE()$value, width=0.2)) +
            geom_errorbar(data = Ft(), aes(x=as.numeric(as.character(Ft()$variable)), y=value, ymax = Ft()$value + FSE()$value, ymin=Ft()$value - FSE()$value, width=0.2)) +
            geom_errorbar(data = Ht(), aes(x=as.numeric(as.character(Ht()$variable)), y=value, ymax = Ht()$value + HSE()$value, ymin=Ht()$value - HSE()$value, width=0.2)) +
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
