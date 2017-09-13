options(shiny.sanitize.errors = FALSE)
library(ggplot2)
library(reshape2)
library(plyr)
library(shinyBS)
library(shinythemes)
library(DT)
library(shiny)
library(emoGG)

#This lets you hit return in place of click
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

#this makes the buttons and input
ui <- tagList(
        tags$head(
          tags$style(HTML("
                                  #tableF {
                                    width:auto;
                                    word-wrap: break-word;
                                    word-break: break-all;
                                    font-family: Courier New;
                                  }
                                  #tableF td {
                                    min-width: 150px;
                                  }
                                  
                                  ")),
          #supress error messages:
          tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
          )
        ),
        bsModal("modal", "Search Results", "", size = "large", DT::dataTableOutput('searchTable')),
        bsModal("converted", "Conversion Results", "", size = "large", DT::dataTableOutput('conversionTable')),
        navbarPage(theme = shinytheme("cerulean"),id = "inTabset", title="StellaTx",
          sidebarPanel(
            tags$head(tags$script(HTML(jscode))),
            h5("Search for a gene"),
            tagAppendAttributes(
              textInput('search', label=NULL, value = '',width = NULL, placeholder = "e.g. 'TCF'"),
               `data-proxy-click` = "go"
             ),
            actionButton("go", "Search!"),
            h4(""),
            tags$head(tags$script(HTML(jscode))),
            tagAppendAttributes(
               textInput('gene1', 'Input NvERTx number', value = "", width = NULL, placeholder =  'e.g. NvERTx.4.100038'),
              `data-proxy-click` = "do"
             ),
            tags$head(tags$script(HTML(jscode))),
            tagAppendAttributes(
              textInput('gene2', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.4.100038'),
              `data-proxy-click` = "do"
            ),
            tags$head(tags$script(HTML(jscode))),
            tagAppendAttributes(
              textInput('gene3', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.4.100038'),
              `data-proxy-click` = "do"
            ),
            tags$head(tags$script(HTML(jscode))),
            tagAppendAttributes(
              textInput('gene4', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.4.100038'),
              `data-proxy-click` = "do"
            ),
            tags$head(tags$script(HTML(jscode))),
            tagAppendAttributes(
            textInput('gene5', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.4.100038'),
              `data-proxy-click` = "do"
            ),
            checkboxInput('log', 'Check for Log2', value = FALSE, width = NULL),
            actionButton("do", "Evaluate!"),
            checkboxInput('returnpdf', 'Output Regen pdf?', FALSE),
            conditionalPanel(
              condition = "input.returnpdf == true",
              strong("PDF size (inches):"),
              sliderInput(inputId="w", label = "width:", min=3, max=20, value=12, width=100, ticks=F),
              sliderInput(inputId="h", label = "height:", min=3, max=20, value=6, width=100, ticks=F),
              br(),
              downloadButton('downloadPlot', 'Download Plot')
            ),
            checkboxInput('returnpdf2', 'Output Embryo pdf?', FALSE),
            conditionalPanel(
              condition = "input.returnpdf2 == true",
              strong("PDF size (inches):"),
              sliderInput(inputId="w", label = "width:", min=3, max=20, value=12, width=100, ticks=F),
              sliderInput(inputId="h", label = "height:", min=3, max=20, value=6, width=100, ticks=F),
              br(),
              downloadButton('downloadPlot2', 'Download Plot')
              ),
            h5("Convert NvERTx.2 to NvERTx.4"),
            tags$head(tags$script(HTML(jscode))),
            tagAppendAttributes(
              textAreaInput('NveConvert', 'Input NvERtx.2 IDs separated by comma:', value = "", placeholder = "e.g: NvERTx.2.133024,NvERTx.2.1000"),
              `data-proxy-click` = "convert"
            ),
            actionButton("convert", "Convert!")
          ),
          tabPanel("Home",
            h2("StellaTx - An embryogenesis & regeneration gene expression plotter"),
            p("Welcome to NvERTx, an embryonic and regenerative transcriptome exploration tool. To learn more about the assembly check out the About Page"),
            p("To get started you can enter a search term (gene name or JGI ID) in the sidebar or if you have a sequence you can BLAST it using the BLAST tool. You can also explore the mFuzz transcript clusters."),
            img(src="EmbryoRegen.png", height = 400, width = 500)
          ),
          tabPanel("Results",
            mainPanel(
              h2("Results:"),
              div(
                tabsetPanel(type = "tabs", 
                  tabPanel("Plots",
                    h4("Regeneration Expression"),
                    plotOutput("plot1"),
                    h4("Embryonic Expression"),
                    plotOutput("plot2")
                  ),
                  tabPanel("Count Tables", 
                    h4("Regeneration average counts (hours post amputation)"),
                    tableOutput('table'),
                    h4("Warner et al. average counts (hours post fertilization)"),
                    tableOutput('table3'),
                    h4("Fischer et al. counts (hours post fertilization)"),
                    tableOutput('table4'),
                    h4("Helm et al. counts (hours post fertilization)"),
                    tableOutput('table5')
                  ),
                  tabPanel("Annotation", 
                    h4("Annotation"),
                    p("Click '+' for extra hits"),
                    DT::dataTableOutput('table2')
                  ),
                  tabPanel("Fasta", 
                    h4("Fasta"),
                    tableOutput("tableF")
                  ),
                  tabPanel("Bibliography", 
                    h4("Literature Hits"),
                    DT::dataTableOutput("tableP")
                  )
                ), class="span7")
              )
            ),
            tabPanel("mFuzz Clusters",
              mainPanel(
                uiOutput("mfuzzModal"),
                tabsetPanel(type = "tabs", 
                  tabPanel("Regeneration Clusters",
                    actionButton("R1", img(src="MfuzzR-1.png", height = 125, width = 125)),
                    actionButton("R2", img(src="MfuzzR-2.png", height = 125, width = 125)),
                    actionButton("R3", img(src="MfuzzR-3.png", height = 125, width = 125)),
                    actionButton("R4", img(src="MfuzzR-4.png", height = 125, width = 125)),
                    actionButton("R5", img(src="MfuzzR-5.png", height = 125, width = 125)),
                    actionButton("R6", img(src="MfuzzR-6.png", height = 125, width = 125)),
                    actionButton("R7", img(src="MfuzzR-7.png", height = 125, width = 125)),
                    actionButton("R8", img(src="MfuzzR-8.png", height = 125, width = 125)),
                    actionButton("R9", img(src="MfuzzR-9.png", height = 125, width = 125))
                  ),
                  tabPanel("Embryonic Clusters",
                    actionButton("E1", img(src="MfuzzE-1.png", height = 125, width = 125)),
                    actionButton("E2", img(src="MfuzzE-2.png", height = 125, width = 125)),
                    actionButton("E3", img(src="MfuzzE-3.png", height = 125, width = 125)),
                    actionButton("E4", img(src="MfuzzE-4.png", height = 125, width = 125)),
                    actionButton("E5", img(src="MfuzzE-5.png", height = 125, width = 125)),
                    actionButton("E6", img(src="MfuzzE-6.png", height = 125, width = 125)),
                    actionButton("E7", img(src="MfuzzE-7.png", height = 125, width = 125)),
                    actionButton("E8", img(src="MfuzzE-8.png", height = 125, width = 125))
                  )
                )
              )
            ),
            tabPanel("Volcano Plots",
              mainPanel(
                selectInput('volcano',label ="Select Comparison", 
                  choices=c("Regeneration_0-2Hpa",
                            "Regeneration_0-4Hpa",
                            "Regeneration_0-8Hpa",
                            "Regeneration_0-12Hpa",
                            "Regeneration_0-16Hpa",
                            "Regeneration_0-20Hpa",
                            "Regeneration_0-24Hpa",
                            "Regeneration_0-36Hpa",
                            "Regeneration_0-48Hpa",
                            "Regeneration_0-60Hpa",
                            "Regeneration_0-72Hpa",
                            "Regeneration_0-96Hpa",
                            "Regeneration_0-120Hpa",
                            "Regeneration_0-144Hpa",
                            "Embryo_Warner_24-48Hpf",
                            "Embryo_Warner_24-72Hpf",
                            "Embryo_Warner_24-96Hpf",
                            "Embryo_Warner_24-120Hpf",
                            "Embryo_Warner_24-144Hpf",
                            "Embryo_Warner_24-168Hpf",
                            "Embryo_Warner_24-192Hpf",
                            "Embryo_Helm_7-12Hpf",
                            "Embryo_Helm_7-24Hpf",
                            "Embryo_Helm_7-120Hpf",
                            "Embryo_Helm_7-240Hpf",
                            "Embryo_Fischer_7-8Hpf",
                            "Embryo_Fischer_7-9Hpf",
                            "Embryo_Fischer_7-10Hpf",
                            "Embryo_Fischer_7-11Hpf",
                            "Embryo_Fischer_7-12Hpf",
                            "Embryo_Fischer_7-13Hpf",
                            "Embryo_Fischer_7-14Hpf",
                            "Embryo_Fischer_7-15Hpf",
                            "Embryo_Fischer_7-16Hpf",
                            "Embryo_Fischer_7-17Hpf",
                            "Embryo_Fischer_7-18Hpf",
                            "Embryo_Fischer_7-19Hpf")
                  ),
                div(style="display:inline-block",
                  sliderInput('cutoff', label="-log(FDR) cutoff",1,100,1,step=0.01, width="200px")
                ),
                div(style="display:inline-block",
                  sliderInput('cpmCut', label="log(CPM) cutoff",0,10,0, width="200px")
                ),
                plotOutput("Vplot", click = "plot_click"
                                      #  dblclick = "plot_dblclick",
                                      # hover = "plot_hover",
                                      #  brush = "plot_brush"
                ),
                h4("Selected Points"),
                dataTableOutput("info")
                           #h4("Extra Info"),
                           #tableOutput("volcSelected")
                )
              ),
              tabPanel("Blast",
                #This block gives us all the inputs:
                mainPanel(
                  p("This only works locally for now...",style="color:red"),
                  p("Click 'BLAST!' to see a sample result.",style="color:red"),
                  textAreaInput('query', 'Input sequence:', value = "", placeholder = "", width = "600px", height="200px"),
                  selectInput("db", "Databse:", choices=c("NvERTx.4","nr"), width="120px"),
                  div(style="display:inline-block",
                    selectInput("program", "Program:", choices=c("blastn","tblastn"), width="100px")),
                  div(style="display:inline-block",
                    selectInput("eval", "e-value:", choices=c(1,0.001,1e-4,1e-5,1e-10), width="120px")),
                  actionButton("blast", "BLAST!")
                ),
                #Basic results output
                mainPanel(
                  h4("Blast Results"),
                  DT::dataTableOutput("blastResults"),
                  p("Alignment:", tableOutput("clicked") ),
                  verbatimTextOutput("alignment")
                )
              ),
              tabPanel("Help",
                mainPanel(
                  h3("How to use this site:"),
                  p("Blah-Blah-Blah-Blah-Blah-Blah-Blah-Blah-")
                )
              ),
              tabPanel("About",
                mainPanel(
                  h3("StellaTx"),
                  p("Understanding the relationship between embryogenesis and regeneration is a long standing question in regenerative biology as both developmental strategies lead to fully functional organisms. Modern sequencing techniques have allowed us to re-examine this question at the systems level. This site provides a set of tools to analyze gene expression in both embryogenesis and regeneration of the starlet sea anemone Nematostella vectensis (Anthozoa, Cnidaria). Here you will find tools to mine RNAseq datasets for embryonic development and oral regeneration of Nematostella vectensis. The data can be accessed by searching for your gene of interest using a gene name, NvERTx.4 ID, or JGI Nemve1 accession number. If you have a sequence you can use the BLAST function to retrieve corresponding NvERTx transcripts. You can also explore a list of genes corresponding to a given embryonic or regeneration co-expression cluster.
                    References for regeneration and novel embryonic datasets: Warner et al(a), 2017, Warner et al(b), 2017. Re-analyzed embryonic datasets: Fischer et al. 2014,Tulin et al. 2013,and Helm et al. 2013"),
                  h3("Data Availability"),
                  p("Datasets are publicly available. The NvERTx.4 trinity assembly is available for download using this direct download link. Raw sequencing reads from Fischer et al. 2014 (0-19Hpf hourly) and Tulin et al. 2013 (0,6,12,18,24Hpf) are available on the Woods Hole Open Access Server here and here. Raw sequencing reads from Helm et al. 2013 (2,7,12,24,120,240 Hpf) are available on the NCBI sequence read archive project number PRJNA189768"),
                  h3("About this site"),
                  p("This site was developed using the Shiny framework and coded by Jake Warner. To report bugs or request data, send us an email at: jwarner\at/unice.fr"),
                  h3("Funding"),
                  p("Funding for this site and it's developers was graciously provided by:"),
                  img(src="arc.png",height=125, width=175),
                  img(src="ircan.png",height=175, width=175),
                  img(src="INSERM.jpeg",height=125, width=175),
                  img(src="UCA.jpeg",height=175, width=175)
                )
              ),
              selected = "Home")
)