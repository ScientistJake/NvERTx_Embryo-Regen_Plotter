options(shiny.sanitize.errors = FALSE)
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
library(RSQLite)
library(DBI)
library(httr)

server <- function(input, output, session){
  con = dbConnect(RSQLite::SQLite(), dbname="db.sqlite3")
  observeEvent(input$do, {
    updateTabsetPanel(session, "inTabset",
                      selected = "Results")
  })

  ##########
  ##########
  ##########
  ##########              ******* Search Function *********
  ##########
  ##########
  ##########
  
  searchTerm <- eventReactive(input$go,{
    trimws(input$search)
  },ignoreNULL= T)
  
  output$searchTable = renderDataTable({
    req(searchTerm())
    search <- dbGetQuery(con, paste0("select * from ER_plotter_annotation where ((Nemve1_tophit LIKE ",paste("'%",as.character(searchTerm()),"%'",sep=""),") OR (Uniprot_ID LIKE ",paste("'%",as.character(searchTerm()),"%'",collapse=", ",sep=""),") OR (Top_nr_hit_eval LIKE ",paste("'%",as.character(searchTerm()),"%'",collapse=", ",sep=""),") OR (Uniprot_Description LIKE ",paste("'%",as.character(searchTerm()),"%'",collapse=", ",sep=""),"))"))
    NvERtx_ID <- search$nvertx_id
    search <- search[c(2,9,10,11)]
    datatable( cbind(' ' = '&oplus;',NvERtx_ID, search),
               rownames = T,
               escape = F,
               options = list(dom = 'ft', 
                              columnDefs = list(
                                list(visible = FALSE, targets = c(0,6)),
                                list(orderable = FALSE, className = 'details-control', targets = 1)
                              )
               ),
               callback = JS("table.column(1).nodes().to$().css({cursor: 'pointer'});
                             var format = function(d) {
                             return '<div style=\"background-color:#eee; padding: .5em;\"> Extra Hits: ' +
                             d[6] + '</div>';
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
    req(searchTerm())
    search <- dbGetQuery(con, paste0("select * from ER_plotter_annotation where ((Nemve1_tophit LIKE ",paste("'%",as.character(searchTerm()),"%'",sep=""),") OR (Uniprot_ID LIKE ",paste("'%",as.character(searchTerm()),"%'",collapse=", ",sep=""),") OR (Top_nr_hit_eval LIKE ",paste("'%",as.character(searchTerm()),"%'",collapse=", ",sep=""),") OR (Uniprot_Description LIKE ",paste("'%",as.character(searchTerm()),"%'",collapse=", ",sep=""),"))"))
    NvERtx_ID <- search$nvertx_id
    search <- search[c(2,9,10,11)]
    row.names(search) <- NvERtx_ID
    Aclicked = input$searchTable_rows_selected
    updateTextInput(session, "gene1", value = row.names(search)[Aclicked[1]])
    updateTextInput(session, "gene2", value = row.names(search)[Aclicked[2]])
    updateTextInput(session, "gene3", value = row.names(search)[Aclicked[3]])
    updateTextInput(session, "gene4", value = row.names(search)[Aclicked[4]])
    updateTextInput(session, "gene5", value = row.names(search)[Aclicked[5]])
  })
  
  observeEvent(input$go, {
    toggleModal(session, "modal", "open")
  })
  
  ##########
  ##########
  ##########
  ##########              ******* MFUZZ Panel *********
  ##########
  ##########
  ##########
  
  ##########
  ##
  ##
  ## Regen Fuzz
  ##
  ##
  ##########
  
  M1 <- eventReactive(input$M1,{
    M1 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-1')")
    row.names(M1) <- M1$nvertx_id
    M1 <- M1[,c(4,5,9,10)]
    M1[order(-M1$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM1 <- renderDataTable({
    req(M1())
    datatable(M1(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M1, {
    toggleModal(session, "m1", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M1())
    m1clicked = input$tableM1_rows_selected
    updateTextInput(session, "gene1", value = row.names(M1())[m1clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M1())[m1clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M1())[m1clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M1())[m1clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M1())[m1clicked[5]])
  })
  
  M2 <- eventReactive(input$M2,{
    M2 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-2')")
    row.names(M2) <- M2$nvertx_id
    M2 <- M2[,c(4,5,9,10)]
    M2[order(-M2$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM2 <- renderDataTable({
    req(M2())
    datatable(M2(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M2, {
    toggleModal(session, "m2", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M2())
    m2clicked = input$tableM2_rows_selected
    updateTextInput(session, "gene1", value = row.names(M2())[m2clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M2())[m2clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M2())[m2clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M2())[m2clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M2())[m2clicked[5]])
  })
  
  M3 <- eventReactive(input$M3,{
    M3 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-3')")
    row.names(M3) <- M3$nvertx_id
    M3 <- M3[,c(4,5,9,10)]
    M3[order(-M3$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM3 <- renderDataTable({
    req(M3())
    datatable(M3(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M3, {
    toggleModal(session, "m3", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M3())
    m3clicked = input$tableM3_rows_selected
    updateTextInput(session, "gene1", value = row.names(M3())[m3clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M3())[m3clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M3())[m3clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M3())[m3clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M3())[m3clicked[5]])
  })
  
  M4 <- eventReactive(input$M4,{
    M4 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-4')")
    row.names(M4) <- M4$nvertx_id
    M4 <- M4[,c(4,5,9,10)]
    M4[order(-M4$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM4 <- renderDataTable({
    req(M4())
    datatable(M4(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M4, {
    toggleModal(session, "m4", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M4())
    m4clicked = input$tableM4_rows_selected
    updateTextInput(session, "gene1", value = row.names(M4())[m4clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M4())[m4clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M4())[m4clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M4())[m4clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M4())[m4clicked[5]])
  })
  
  M5 <- eventReactive(input$M5,{
    M5 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-5')")
    row.names(M5) <- M5$nvertx_id
    M5 <- M5[,c(4,5,9,10)]
    M5[order(-M5$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM5 <- renderDataTable({
    req(M5())
    datatable(M5(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M5, {
    toggleModal(session, "m5", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M5())
    m5clicked = input$tableM5_rows_selected
    updateTextInput(session, "gene1", value = row.names(M5())[m5clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M5())[m5clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M5())[m5clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M5())[m5clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M5())[m5clicked[5]])
  })
  
  M6 <- eventReactive(input$M6,{
    M6 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-6')")
    row.names(M6) <- M6$nvertx_id
    M6 <- M6[,c(4,5,9,10)]
    M6[order(-M6$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM6 <- renderDataTable({
    req(M6())
    datatable(M6(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M6, {
    toggleModal(session, "m6", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M6())
    m6clicked = input$tableM6_rows_selected
    updateTextInput(session, "gene1", value = row.names(M6())[m6clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M6())[m6clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M6())[m6clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M6())[m6clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M6())[m6clicked[5]])
  })
  
  M7 <- eventReactive(input$M7,{
    M7 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-7')")
    row.names(M7) <- M7$nvertx_id
    M7 <- M7[,c(4,5,9,10)]
    M7[order(-M7$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM7 <- renderDataTable({
    req(M7())
    datatable(M7(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M7, {
    toggleModal(session, "m7", "open")
  })
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M7())
    m7clicked = input$tableM7_rows_selected
    updateTextInput(session, "gene1", value = row.names(M7())[m7clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M7())[m7clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M7())[m7clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M7())[m7clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M7())[m7clicked[5]])
  })
  
  M8 <- eventReactive(input$M8,{
    M8 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-8')")
    row.names(M8) <- M8$nvertx_id
    M8 <- M8[,c(4,5,9,10)]
    M8[order(-M8$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM8 <- renderDataTable({
    req(M8())
    datatable(M8(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M8, {
    toggleModal(session, "m8", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M8())
    m8clicked = input$tableM8_rows_selected
    updateTextInput(session, "gene1", value = row.names(M8())[m8clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M8())[m8clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M8())[m8clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M8())[m8clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M8())[m8clicked[5]])
  })
  
  M9 <- eventReactive(input$M9,{
    M9 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_R_Clust in ('R-9')")
    row.names(M9) <- M9$nvertx_id
    M9 <- M9[,c(4,5,9,10)]
    M9[order(-M9$Mfuzz_R_Score),]
  },ignoreNULL= T)
  
  output$tableM9 <- renderDataTable({
    req(M9())
    datatable(M9(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$M9, {
    toggleModal(session, "m9", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(M9())
    m9clicked = input$tableM9_rows_selected
    updateTextInput(session, "gene1", value = row.names(M9())[m9clicked[1]])
    updateTextInput(session, "gene2", value = row.names(M9())[m9clicked[2]])
    updateTextInput(session, "gene3", value = row.names(M9())[m9clicked[3]])
    updateTextInput(session, "gene4", value = row.names(M9())[m9clicked[4]])
    updateTextInput(session, "gene5", value = row.names(M9())[m9clicked[5]])
  })
  
  #  M10 <- eventReactive(input$M10,{
  #    M10 <- annotations[annotations$Mfuzz_R_Clust %in% c("R-10"),c(4,5,9,10)]
  #    M10[order(-M10$Mfuzz_R_Score),]
  #  },ignoreNULL= T)
  #  
  #  output$tableM10 <- renderDataTable({
  #    M10()
  #  })
  #  
  #  observeEvent(input$M10, {
  #    toggleModal(session, "m10", "open")
  #  })
  #  
  #  #this code block updates the inputs if a row is clicked on the pop-up table
  #  observe({
  #    m10clicked = input$tableM10_rows_selected
  #    updateTextInput(session, "gene1", value = row.names(M10())[m10clicked[1]])
  #    updateTextInput(session, "gene2", value = row.names(M10())[m10clicked[2]])
  #    updateTextInput(session, "gene3", value = row.names(M10())[m10clicked[3]])
  #    updateTextInput(session, "gene4", value = row.names(M10())[m10clicked[4]])
  #    updateTextInput(session, "gene5", value = row.names(M10())[m10clicked[5]])
  #  })
  
  ##########
  ##
  ##
  ## Embryo Fuzz
  ##
  ##
  ##########
  
  E1 <- eventReactive(input$E1,{
    E1 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-1')")
    row.names(E1) <- E1$nvertx_id
    E1 <- E1[,c(6,7,9,10)]
    E1[order(-E1$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE1 <- renderDataTable({
    req(E1())
    datatable(E1(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E1, {
    toggleModal(session, "e1", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E1())
    e1clicked = input$tableE1_rows_selected
    updateTextInput(session, "gene1", value = row.names(E1())[e1clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E1())[e1clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E1())[e1clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E1())[e1clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E1())[e1clicked[5]])
  })
  
  E2 <- eventReactive(input$E2,{
    E2 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-2')")
    row.names(E2) <- E2$nvertx_id
    E2 <- E2[,c(6,7,9,10)]
    E2[order(-E2$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE2 <- renderDataTable({
    req(E2())
    datatable(E2(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E2, {
    toggleModal(session, "e2", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E2())
    e2clicked = input$tableE2_rows_selected
    updateTextInput(session, "gene1", value = row.names(E2())[e2clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E2())[e2clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E2())[e2clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E2())[e2clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E2())[e2clicked[5]])
  })
  
  E3 <- eventReactive(input$E3,{
    E3 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-3')")
    row.names(E3) <- E3$nvertx_id
    E3 <- E3[,c(6,7,9,10)]
    E3[order(-E3$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE3 <- renderDataTable({
    req(E3())
    datatable(E3(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E3, {
    toggleModal(session, "e3", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E3())
    e3clicked = input$tableE3_rows_selected
    updateTextInput(session, "gene1", value = row.names(E3())[e3clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E3())[e3clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E3())[e3clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E3())[e3clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E3())[e3clicked[5]])
  })
  
  E4 <- eventReactive(input$E4,{
    E4 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-4')")
    row.names(E4) <- E4$nvertx_id
    E4 <- E4[,c(6,7,9,10)]
    E4[order(-E4$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE4 <- renderDataTable({
    req(E4())
    datatable(E4(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E4, {
    toggleModal(session, "e4", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E4())
    e4clicked = input$tableE4_rows_selected
    updateTextInput(session, "gene1", value = row.names(E4())[e4clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E4())[e4clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E4())[e4clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E4())[e4clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E4())[e4clicked[5]])
  })
  
  E5 <- eventReactive(input$E5,{
    E5 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-5')")
    row.names(E5) <- E5$nvertx_id
    E5 <- E5[,c(6,7,9,10)]
    E5[order(-E5$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE5 <- renderDataTable({
    req(E5())
    datatable(E5(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E5, {
    toggleModal(session, "e5", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E5())
    e5clicked = input$tableE5_rows_selected
    updateTextInput(session, "gene1", value = row.names(E5())[e5clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E5())[e5clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E5())[e5clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E5())[e5clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E5())[e5clicked[5]])
  })
  
  E6 <- eventReactive(input$E6,{
    E6 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-6')")
    row.names(E6) <- E6$nvertx_id
    E6 <- E6[,c(6,7,9,10)]
    E6[order(-E6$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE6 <- renderDataTable({
    req(E6())
    datatable(E6(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E6, {
    toggleModal(session, "e6", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E6())
    e6clicked = input$tableE6_rows_selected
    updateTextInput(session, "gene1", value = row.names(E6())[e6clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E6())[e6clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E6())[e6clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E6())[e6clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E6())[e6clicked[5]])
  })
  
  E7 <- eventReactive(input$E7,{
    E7 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-7')")
    row.names(E7) <- E7$nvertx_id
    E7 <- E7[,c(6,7,9,10)]
    E7[order(-E7$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE7 <- renderDataTable({
    req(E7())
    datatable(E7(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E7, {
    toggleModal(session, "e7", "open")
  })
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E7())
    e7clicked = input$tableE7_rows_selected
    updateTextInput(session, "gene1", value = row.names(E7())[e7clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E7())[e7clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E7())[e7clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E7())[e7clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E7())[e7clicked[5]])
  })
  
  E8 <- eventReactive(input$E8,{
    E8 <- dbGetQuery(con, "select * from ER_plotter_annotation where Mfuzz_E_Clust in ('E-8')")
    row.names(E8) <- E8$nvertx_id
    E8 <- E8[,c(6,7,9,10)]
    E8[order(-E8$Mfuzz_E_Score),]
  },ignoreNULL= T)
  
  output$tableE8 <- renderDataTable({
    req(E8())
    datatable(E8(),extensions = 'Buttons', options = list(dom = 'Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })
  
  observeEvent(input$E8, {
    toggleModal(session, "e8", "open")
  })
  
  #this code block updates the inputs if a row is clicked on the pop-up table
  observe({
    req(E8())
    e8clicked = input$tableE8_rows_selected
    updateTextInput(session, "gene1", value = row.names(E8())[e8clicked[1]])
    updateTextInput(session, "gene2", value = row.names(E8())[e8clicked[2]])
    updateTextInput(session, "gene3", value = row.names(E8())[e8clicked[3]])
    updateTextInput(session, "gene4", value = row.names(E8())[e8clicked[4]])
    updateTextInput(session, "gene5", value = row.names(E8())[e8clicked[5]])
  })

  
  ##########
  ##########
  ##########
  ##              ******* Results Panel *********
  ##########
  ##########
  ##########
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['gene1']])) {
      updateTextInput(session, "gene1", value = query[['gene1']])
    }
  })
 
  #this assigns the NvERTx numbers to the nve object
  nve <- eventReactive(input$do, {
    nve <- c(trimws(input$gene1), trimws(input$gene2), trimws(input$gene3), trimws(input$gene4), trimws(input$gene5))
    nve[nve != ""]
  }, ignoreNULL= T)
  
  ##########
  ##
  ##
  ## Regen Tables
  ##
  ##
  ##########
  
  #this subsets the count table by the nve numbers
  gene <- reactive({
    req(nve())
    if (as.numeric(input$log) >0 ){
      counts <- dbGetQuery(con, paste0("select * from ER_plotter_regen_cpm where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
      rownames(counts) <- counts$nvertx_id
      ID <- counts$nvertx_id
      counts <- counts[-c(1)]
      counts <- log(counts +1 , base=2)
      counts$ID <- ID
      names(counts) <- c("-1","0","2","4","8","12","16","20","24","36","48","60","72","96","120","144","ID")
      counts<- counts[complete.cases(counts),]
    } else {
      counts <- dbGetQuery(con, paste0("select * from ER_plotter_regen_cpm where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
      rownames(counts) <- counts$nvertx_id
      ID <- counts$nvertx_id
      counts <- counts[-c(1)]
      counts$ID <- ID
      names(counts) <- c("-1","0","2","4","8","12","16","20","24","36","48","60","72","96","120","144","ID")
      counts<- counts[complete.cases(counts),]
    }
    counts
  })

  #same workflow, subsetting the SE values
  geneSE <- reactive({
    req(nve())
    if (as.numeric(input$log) >0 ){
      req(nve())
      countsSE <- dbGetQuery(con, paste0("select * from ER_plotter_regen_log_se where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
      rownames(countsSE) <- countsSE$nvertx_id
      ID <- countsSE$nvertx_id
      countsSE <- countsSE[-c(1)]
      countsSE$ID <- ID
      names(countsSE) <- c("-1","0","2","4","8","12","16","20","24","36","48","60","72","96","120","144","ID")
      countsSE<- countsSE[complete.cases(countsSE),]
    } else {
      countsSE <- dbGetQuery(con, paste0("select * from ER_plotter_regen_se where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
      rownames(countsSE) <- countsSE$nvertx_id
      ID <- countsSE$nvertx_id
      countsSE <- countsSE[-c(1)]
      countsSE$ID <- ID
      names(countsSE) <- c("-1","0","2","4","8","12","16","20","24","36","48","60","72","96","120","144","ID")
      countsSE<- countsSE[complete.cases(countsSE),]
    }
    countsSE
  })
  
  output$table <- renderTable({
    req(gene())
    gene()[,c(ncol(gene()),1:ncol(gene())-1)]
  })
  
  
  ##########
  ##
  ##
  ## Regen Plots
  ##
  ##
  ##########
  
  #melt the data for plotting
  geneMolten <- reactive({
    req(gene())
    melt(gene())
  }) #this puts samples as rows, genes as columns 
  geneMoltenSE <- reactive({
    req(geneSE())
    melt(geneSE())
  }) #this puts samples as rows, genes as columns 
  
  # this uses the SE values to make the error bar limits
  limits <- reactive({
    req(geneMolten(),geneMoltenSE())
    limits <- aes(ymax = geneMolten()$value + geneMoltenSE()$value, ymin=geneMolten()$value - geneMoltenSE()$value)
    limits
    }) # This is the calculation for the error bars
  #for whatever reason I couldnt do the same for the embryonic. I had to explicitly map them in the Geom_errorbar() below
  
  p <- reactive({
    req(geneMolten())
    ggplot(geneMolten(), aes(x=as.numeric(as.character(geneMolten()$variable)), y=value, colour=ID)) + geom_line() +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
    
  })
  r <- reactive({
    req(geneMolten())
    ggplot(geneMolten(), aes(x=as.numeric(as.character(geneMolten()$variable)), y=value, colour=ID)) + geom_line() +
      geom_emoji(emoji="1f433") +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
  })
  s <- reactive({
    req(geneMolten())
    ggplot(geneMolten(), aes(x=as.numeric(as.character(geneMolten()$variable)), y=value, colour=ID)) + geom_line() +
      geom_emoji(emoji="1f412") +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
  })
  t <- reactive({
    req(geneMolten())
    ggplot(geneMolten(), aes(x=as.numeric(as.character(geneMolten()$variable)), y=value, colour=ID)) + geom_line() +
      geom_emoji(emoji="1f369") +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
  })
  
  # 3% chance to get an emoji plot!
  rplot <- reactive({
    req(r(),s(),t(),p())
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
  
  # change the axis label for the plot if log is ticked
  y_label <- reactive({
    if (as.numeric(input$log) >0 ){
      return(c("Log2(Counts per million +1)"))
    } else {
      return(c("Counts per million"))
    }
  })
  
  output$plot1 <- renderPlot({ 
    req(rplot(),limits(),y_label())
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
      ggsave(file, plot = p()+ scale_x_continuous(minor_breaks = NULL, breaks=c(0,6,12,24,48,72,96,120,144,168,192,120,240)) +
               geom_errorbar(limits(), width=0.2) +
               ylab(y_label()) +
               ylab("Log2(Counts per million +1)") +
               xlab("Hours Post Fertilization"), device = device)
    }
  )
  
  ##########
  ##
  ##
  ## Embryo Tables
  ##
  ##
  ##########
  
  Egene <- reactive({
    req(nve())
    EmbryoSQL <- dbGetQuery(con, paste0("select * from ER_plotter_embryo_cpm where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
    rownames(EmbryoSQL) <- EmbryoSQL$nvertx_id
    EmbryoSQL$ID <- EmbryoSQL$nvertx_id
    EmbryoSQL <- EmbryoSQL[-c(1)]
    names(EmbryoSQL) <- c("24","48","72","96","120","144","168","192","0","1","2","3","4","5","6","7","8","9","10","11",
                          "12","13","14","15","16","17","18","19","2","7","12","24","120","240","2","7","12","24","120","ID")
    EmbryoSQL <- EmbryoSQL[complete.cases(EmbryoSQL),]
    EmbryoSQL
  })

  #this gets the mean
  Egenex <- reactive({
    req(Egene())
    Egene()[c(2:4, 6:10, 12:15,17:20, 22:28,34:40)]
  })
  
  #same workflow, subsetting the SE values
  EgeneSE <- reactive({
    req(nve())
    EmbryoSESQL <- dbGetQuery(con, paste0("select * from ER_plotter_embryo_se where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
    rownames(EmbryoSESQL) <- EmbryoSESQL$nvertx_id
    EmbryoSESQL$ID <- EmbryoSESQL$nvertx_id
    EmbryoSESQL <- EmbryoSESQL[-c(1)]
    names(EmbryoSESQL) <- c("24","48","72","96","120","144","168","192","0","1","2","3","4","5","6","7","8","9","10","11",
                            "12","13","14","15","16","17","18","19","2","7","12","24","120","240","ID")
    #make sure the SEs are the same as the counts
    EmbryoSESQL[which(EmbryoSESQL$ID %in% Egene()$ID),]
  })
  
  Egenet <- reactive({
    req(Egenex())
    melt(Egenex())
  })
  
  #this separates the Embryo datasets so they plot with different shapes
  E <- reactive({
    req(Egene())
    Egene()[c(1:8,40)]
  })
  Et <- reactive({
    req(E())
    melt(E())
  })
  ERSE <- reactive({
    req(EgeneSE())
    ERSEx <- EgeneSE()[c(1:8,35)]
    melt(ERSEx)
  })
  
  Fisch <- reactive({
    req(Egene())
    Egene()[c(9:28,40)]
  })
  Ft <- reactive({
    req(Fisch())
    melt(Fisch())
  })
  FSE <- reactive({
    req(EgeneSE())
    FSEx <- EgeneSE()[c(9:28,35)]
    melt(FSEx)
  })
  
  H <- reactive({
    req(Egene())
    Egene()[c(29:34,40)]
  })
  Ht <- reactive({
    req(H())
    melt(H())
  })
  HSE <- reactive({
    req(EgeneSE())
    HSEx <- EgeneSE()[c(29:34,35)]
    melt(HSEx)
  })
  
  #first we output the table.  It looks nicer in the long format, hence the 't' for transpose
  #tables are reordered with Id first when printed : last column first then column from 1 to last-1
  output$table3 <- renderTable({
    req(E())
    E()[,c(ncol(E()),1:ncol(E())-1)]
  })
  output$table4 <- renderTable({
    req(Fisch())
    Fisch()[,c(ncol(Fisch()),1:ncol(Fisch())-1)]
  })
  output$table5 <- renderTable({
    req(H())
    H()[,c(ncol(H()),1:ncol(H())-1)]
  })
  
  ##########
  ##
  ##
  ## Embryo Plots
  ##
  ##
  ##########
  
  q <- reactive({
    req(Egenet())
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
    req(q())
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
  
  ##########
  ##
  ##
  ## annotations
  ##
  ##
  ##########
  
  annot <- reactive({
    req(nve())
    dbGetQuery(con, paste0("select * from ER_plotter_annotation where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
  })

  # build up the annotation table with clickable links:
  annotTable <- reactive({
    req(annot())
    uprot <- annot()$Uniprot_ID
    uprotURL <- lapply(uprot,function(i) {
      if (is.na(i)) {
        print(c("Not Available"))
      } else {
        if (i=="No_Uniprotmatch") {
          print(c("No Uniprot match"))
        } else {
          print(paste0("<a href='http://www.uniprot.org/uniprot/",i,"' target='_blank'>",i,"</a>"))
        }}})
    
    annotTable <- annot()
    annotTable$Uniprot_ID <- uprotURL
    #
    NCBI <- strsplit(as.character(annot()$Top_nr_hit_eval), '|', fixed=TRUE)
    NCBI1 <- lapply(NCBI, '[', 1)
    NCBI2 <- lapply(NCBI, '[', 2)
    NCBI3 <- lapply(NCBI, '[', 3)
    NCBI4 <- lapply(NCBI, '[', 4)
    NCBI5 <- lapply(NCBI, '[', 5)
    NCBILink <- lapply(NCBI4,function(i) {
      if (is.na(i)) {
        print(c("No_Nr_Hit"))
      } else {
        print(paste0("<a href='https://www.ncbi.nlm.nih.gov/protein/",i,"' target='_blank'>",i,"</a>"))
      }})
    
    annotTable$Top_nr_hit_eval <- paste0(NCBI1,NCBI2,NCBI3,NCBILink,NCBI5)
    annotTable
  })
  
  #this makes the datatable expandle when you click the icon
  output$table2 <- DT::renderDataTable({
    req(annotTable())
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
  
  ##########
  ##
  ##
  ## Pub Med Hits
  ##
  ##
  ##########
  
  ##output the PubMed URLs:
  pubURLs <- reactive({
    req(annot())
    geneIDs <- strsplit(as.character(annot()$Top_nr_hit_eval), '|', fixed=TRUE)
    uid <- lapply(geneIDs, '[', 2)
    pid <- lapply(geneIDs, '[', 4)
    
    #This gets the primary reference
    urls <- lapply(uid,function(i) {
      if (is.na(i)) {
        print(c("Not available"))
      } else {
        print(paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed?LinkName=protein_pubmed&from_uid=",i,"' target='_blank'>Primary PubMed Link</a>"))
      }})
    
    #This tries to get the link to more references:
    
    more <- lapply(uid,function(i) {
      if (is.na(i)) {
        print(c("Not available"))
      } else {
        pubExtra <- getURL(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&dbto=pubmed&cmd=neighbor&retmode=xml&id=",i))
        pubExtraParse <- xmlParse(pubExtra)
        xml_data <- xmlToList(pubExtraParse)
        pubExtraID <- as.list(xml_data[["LinkSet"]][["LinkSetDb"]][["Link"]][["Id"]])
        if (length(pubExtraID) > 0){
          print(paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed?linkname=pubmed_pubmed&from_uid=",pubExtraID,"' target='_blank'>More PubMed Links</a>"))
        } else {
          print(c("Not available"))
        }
      }})
    
    pblast <- lapply(pid,function(i) {
      if (is.na(i)) {
        print(c("Not available"))
      } else {
        print(paste0("<a href='http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=",i," 'target='_blank'>Paper Blast Link</a>"))
      }})
    
    urltable <- cbind(annot()$nvertx_id, urls, more,pblast)
    colnames(urltable) <- c("NvERTx_ID", "Primary PubMed Link", "Similar PubMed Links", "PaperBlast")
    urltable
  })
  
  output$tableP <- DT::renderDataTable({
    req(pubURLs())
    datatable(pubURLs(), 
              options = list(dom = 't'))
  })
  
  ##########
  ##
  ##
  ## Fasta
  ##
  ##
  ##########
  
  output$tableF <- renderTable({
    req(nve())
    fastaSQL <- dbGetQuery(con, paste0("select * from ER_plotter_fasta where nvertx_id in (",paste("'",as.character(nve()),"'",collapse=", ",sep=""),")"))
    rownames(fastaSQL) <- fastaSQL$nvertx_id
    fastaSQL <- fastaSQL[-c(1)]
    names(fastaSQL) <- c("Sequence")
    fastaSQL
  },rownames = TRUE)
  
  
  ##########
  ##########
  ##########
  ##########              ******* Conversion tool *********
  ##########
  ##########
  ##########
  
  IDtable <- eventReactive(input$convert == 1,{
      IDs <- read.table(file="IDtable.txt", sep='\t', header=T, quote="")
      IDs
  },ignoreNULL= T)
  
  NewIDS <- eventReactive(input$convert,{
    old_IDS <- trimws(unlist(strsplit(as.character(input$NveConvert), ',', fixed=TRUE)))
    IDtable()[IDtable()$NvERTx.2 %in% old_IDS,]
  },ignoreNULL= T)
  
  output$conversionTable <- renderDataTable({
    req(NewIDS())
    NewIDS()
  })
  
  observeEvent(input$convert, {
    toggleModal(session, "converted", "open")
  })
  
  observe({
    req(NewIDS())
    NewIDSclicked = input$conversionTable_rows_selected
    updateTextInput(session, "gene1", value = NewIDS()$NvERTx.4[NewIDSclicked[1]])
    updateTextInput(session, "gene2", value = NewIDS()$NvERTx.4[NewIDSclicked[2]])
    updateTextInput(session, "gene3", value = NewIDS()$NvERTx.4[NewIDSclicked[3]])
    updateTextInput(session, "gene4", value = NewIDS()$NvERTx.4[NewIDSclicked[4]])
    updateTextInput(session, "gene5", value = NewIDS()$NvERTx.4[NewIDSclicked[5]])
  })
  
  ##########
  ##########
  ##########
  ##              ******* BLAST Panel *********
  ##########
  ##########
  ##########
  
  
  blastresults <- eventReactive(input$blast, {
    
    #gather input and set up temp file
    query <- input$query
    tmp <- tempfile(fileext = ".fa")
    
    #if else chooses the right database
    if (input$db =="NvERTx.4"){
      db <- c("/Users/Jake/Documents/GitHub_repos/NvER_plotter_django/nemVec_ER/blast_db/NvERTx.4")
      remote <- c("")
    } else {
      db <- c("nr")
      #add remote option for nr since we don't have a local copy
      remote <- c("-remote")
    }
    
    #this makes sure the fasta is formatted properly
    if (startsWith(query, ">")){
      writeLines(query, tmp)
    } else {
      writeLines(paste0(">Query\n",query), tmp)
    }
    
    #calls the blast
    system(paste0(input$program," -query ",tmp," -db ",db," -evalue ",input$eval," -outfmt 5 -max_hsps 1 -max_target_seqs 10 ",remote), intern = T)
  }, ignoreNULL= T)
  
  #Now to parse the results...
  parsedresults <- reactive({
    req(blastresults())
    xmlparsed <- xmlParse(blastresults())
    xmltop = xmlRoot(xmlparsed)
      
    #the first chunk is for multi-fastas
    results <- xpathApply(xmlparsed, '//Iteration',function(row){
      query_ID <- getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
      hit_IDs <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
      hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
      bitscore <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score') %>% sapply(., xmlValue)
      eval <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue') %>% sapply(., xmlValue)
      cbind(query_ID,hit_IDs,hit_length,bitscore,eval)
    })
    #this ensures that NAs get added for no hits
    results <-  rbind.fill(lapply(results,function(y){as.data.frame((y),stringsAsFactors=FALSE)}))
  })
  
  #makes the datatable
  output$blastResults <- renderDataTable({
    req(parsedresults())
    parsedresults()
  })
  
  #this chunk gets the alignemnt information from a clicked row
  output$clicked <- renderTable({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmlparsed <- xmlParse(blastresults())
      blastclicked = input$blastResults_rows_selected[length(input$blastResults_rows_selected)]
      tableout<- data.frame(parsedresults()[blastclicked,])
      print(tableout)
      names(tableout) <- c("Query ID","Hit ID", "Length", "Bit Score", "e-value")
      data.frame(tableout)
    }
  },rownames =F,colnames =T)
  
  #this chunk makes the alignments for clicked rows
  output$alignment <- renderText({
    if(is.null(input$blastResults_rows_selected)){}
    else{
      xmlparsed <- xmlParse(blastresults())
      xmltop = xmlRoot(xmlparsed)
      
      
      clicked = input$blastResults_rows_selected[length(input$blastResults_rows_selected)]
      
      #loop over the xml to get the alignments
      align <- xpathApply(xmlparsed, '//Iteration',function(row){
        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
        rbind(top,mid,bottom)
      })
      
      #split the alignments every 40 carachters to get a "wrapped look"
      alignx <- do.call("cbind", align)
      splits <- strsplit(gsub("(.{40})", "\\1,", alignx[1:3,clicked]),",")
      
      #paste them together with returns '\n' on the breaks
      split_out <- lapply(1:length(splits[[1]]),function(i){
        rbind(paste0("Q-",splits[[1]][i],"\n"),paste0("M-",splits[[2]][i],"\n"),paste0("H-",splits[[3]][i],"\n"))
      })
      unlist(split_out)
    }
  })
  observe({
    req(parsedresults())
    if (input$db =="NvERTx.4"){
      blastclicked = input$blastResults_rows_selected
      updateTextInput(session, "gene1", value = parsedresults()$hit_IDs[blastclicked[1]])
      updateTextInput(session, "gene2", value = parsedresults()$hit_IDs[blastclicked[2]])
      updateTextInput(session, "gene3", value = parsedresults()$hit_IDs[blastclicked[3]])
      updateTextInput(session, "gene4", value = parsedresults()$hit_IDs[blastclicked[4]])
      updateTextInput(session, "gene5", value = parsedresults()$hit_IDs[blastclicked[5]])
    }
  })
  
  ##########
  ##########
  ##########
  ##              ******* Volcano *********
  ##########
  ##########
  ##########
  
  
  gene_list <- reactive({
    list <- paste0(input$volcano,".txt")
    gene_list <- read.table(file= list, header=T, quote = "", sep='\t')
    gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$FDR < 0.05)
    #gene_listcut <- gene_list[gene_list$FDR < input$cutoff,]
    gene_listcut <- gene_list[gene_list$logCPM > input$cpmCut,]
    gene_listcut$FDR = -log10(gene_listcut$FDR)
    gene_listcut <- gene_listcut[gene_listcut$FDR > input$cutoff,]
    gene_listcut
  })
  ##Construct the plot object
  
  vplot <- reactive({
    req(gene_list())
    ggplot(data=gene_list(), aes(x=logFC, y=FDR, colour=threshold)) +
      geom_point(alpha=0.4, size=1.75) +
      xlim(c(-10, 10)) + ylim(c(0, 65)) +
      theme(axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
            axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
            legend.position = "none"
            ) +
      coord_cartesian() +
      ylab("-log10 FDR") +
      xlab("log2 fold change")
  })
  
  output$Vplot <- renderPlot({ 
    req(vplot())
    #to get the graph to show up in shiny you need to print 
    vplot()
  })
  
  output$info <- DT::renderDataTable({
    #selected_points <- gene_list()[0, ]
    #str(selected_points)
    #selected_points <<- rbind(selected_points, nearPoints(gene_list(), input$plot_click))
    out <- nearPoints(gene_list(), input$plot_click, xvar = "logFC", yvar = "FDR")
    # With base graphics, need to tell it what the x and y variables are.
    get <- rownames(out)
    extras <- dbGetQuery(con, paste0("select nvertx_id, Uniprot_Description, Top_nr_hit_eval from ER_plotter_annotation where nvertx_id in (",paste("'",as.character(get),"'",collapse=", ",sep=""),")"))
    out <- merge(out,extras, by.x="row.names", by.y="nvertx_id")
    row.names(out) <- out$Row.names
    out <- out[c(2,3,4,5,6,8,9)]
    datatable(
      out, extensions = 'Buttons', options = list(dom = 'Btp',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))
    ) %>% formatRound('logFC', 2) %>% formatRound('logCPM', 2) %>% formatRound('FDR', 2)
    # nearPoints() also works with hover and dblclick events
  })
  
  #output$volcSelected <- renderTable({
  #  hits <- nearPoints(gene_list(), input$plot_click, xvar = "logFC", yvar = "FDR")
  #  volcanoclicked = input$info_rows_selected
  #  hitlist <- rownames(hits)[volcanoclicked]
  #  x <- dbGetQuery(con, paste0("select nvertx_id, Nemve1_tophit, Uniprot_Description, Top_nr_hit_eval from ER_plotter_annotation where nvertx_id in (",paste("'",as.character(hitlist),"'",collapse=", ",sep=""),")"))
  #  colnames(x) <- c("NvERTxID", "Nemve1", "Uniprot hit", "Top nr hit")
  #  x
  #})
  
  observe({
    req(vplot())
    hits <- nearPoints(gene_list(), input$plot_click, xvar = "logFC", yvar = "FDR")
      volcanoclicked = input$info_rows_selected
      updateTextInput(session, "gene1", value = rownames(hits)[volcanoclicked[1]])
      updateTextInput(session, "gene2", value = rownames(hits)[volcanoclicked[2]])
      updateTextInput(session, "gene3", value = rownames(hits)[volcanoclicked[3]])
      updateTextInput(session, "gene4", value = rownames(hits)[volcanoclicked[4]])
      updateTextInput(session, "gene5", value = rownames(hits)[volcanoclicked[5]])
  })
  
  #output$info <- renderText({
  # xy_str <- function(e) {
  #    if(is.null(e)) return("NULL\n")
  #    paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
  #  }
  #  xy_range_str <- function(e) {
  #    if(is.null(e)) return("NULL\n")
  #    paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
  #         " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
  #  }
  #  
  #  paste0(
  #    "click: ", xy_str(input$plot_click),
  #    "dblclick: ", xy_str(input$plot_dblclick),
  #    "hover: ", xy_str(input$plot_hover),
  #    "brush: ", xy_range_str(input$plot_brush)
  #  )
  #})
  
  #g = ggplot(data=gene_list, aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  #  geom_point(alpha=0.4, size=1.75) +
    #opts(legend.position = "none") +
  #  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  #  xlab("log2 fold change") + ylab("-log10 p-value")
  #g
  
  session$onSessionEnded(function() {
    dbDisconnect(con)
  })
}
