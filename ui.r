options(shiny.sanitize.errors = FALSE)
library(ggplot2)
library(reshape2)
library(plyr)
library(shinyBS)
library(shinythemes)
library(DT)
library(shiny)
library(emoGG)

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
ui <- fluidPage(theme = shinytheme("cerulean"),
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
                headerPanel('NvERTx Embryo-Regen Plotter'),
                h2( "Now running on SQLite! =) "),
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
                
                mainPanel(
                  bsModal("modal", "Search Results", "", size = "large", DT::dataTableOutput('searchTable')),
                  bsModal("converted", "Conversion Results", "", size = "large", DT::dataTableOutput('conversionTable'))
                ),
                mainPanel(
                  tabsetPanel(type = "tabs", 
                              tabPanel("Regeneration Clusters",
                                       actionButton("M1", img(src="MfuzzR-1.png", height = 125, width = 125)),
                                       bsModal("m1", "Mfuzz1", "", size = "large", img(src="MfuzzR-1.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c1-plot.png", height = 450, width = 600),dataTableOutput('tableM1')),
                                       actionButton("M2", img(src="MfuzzR-2.png", height = 125, width = 125)),
                                       bsModal("m2", "Mfuzz2", "", size = "large", img(src="MfuzzR-2.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c2-plot.png", height = 450, width = 600),dataTableOutput('tableM2')),
                                       actionButton("M3", img(src="MfuzzR-3.png", height = 125, width = 125)),
                                       bsModal("m3", "Mfuzz3", "", size = "large", img(src="MfuzzR-3.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c3-plot.png", height = 450, width = 600),dataTableOutput('tableM3')),
                                       actionButton("M4", img(src="MfuzzR-4.png", height = 125, width = 125)),
                                       bsModal("m4", "Mfuzz4", "", size = "large", img(src="MfuzzR-4.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c4-plot.png", height = 450, width = 600),dataTableOutput('tableM4')),
                                       actionButton("M5", img(src="MfuzzR-5.png", height = 125, width = 125)),
                                       bsModal("m5", "Mfuzz5", "", size = "large", img(src="MfuzzR-5.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c5-plot.png", height = 450, width = 600),dataTableOutput('tableM5')),
                                       actionButton("M6", img(src="MfuzzR-6.png", height = 125, width = 125)),
                                       bsModal("m6", "Mfuzz6", "", size = "large", img(src="MfuzzR-6.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c6-plot.png", height = 450, width = 600),dataTableOutput('tableM6')),
                                       actionButton("M7", img(src="MfuzzR-7.png", height = 125, width = 125)),
                                       bsModal("m7", "Mfuzz7", "", size = "large", img(src="MfuzzR-7.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c7-plot.png", height = 450, width = 600),dataTableOutput('tableM7')),
                                       actionButton("M8", img(src="MfuzzR-8.png", height = 125, width = 125)),
                                       bsModal("m8", "Mfuzz8", "", size = "large", img(src="MfuzzR-8.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c8-plot.png", height = 450, width = 600),dataTableOutput('tableM8')),
                                       actionButton("M9", img(src="MfuzzR-9.png", height = 125, width = 125)),
                                       bsModal("m9", "Mfuzz9", "", size = "large", img(src="MfuzzR-9.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="BP-R-c9-plot.png", height = 450, width = 600),dataTableOutput('tableM9'))
                                       #actionButton("M10", img(src="MfuzzR-10.png", height = 125, width = 125)),
                                       #bsModal("m10", "Mfuzz10", "", size = "large", img(src="MfuzzR-10.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="10-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM10'))
                                      ),
                              tabPanel("Embryonic Clusters",
                                       actionButton("E1", img(src="MfuzzE-1.png", height = 125, width = 125)),
                                       bsModal("e1", "Mfuzz1", "", size = "large", img(src="MfuzzE-1.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="1-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE1')),
                                       actionButton("E2", img(src="MfuzzE-2.png", height = 125, width = 125)),
                                       bsModal("e2", "Mfuzz2", "", size = "large", img(src="MfuzzE-2.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="2-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE2')),
                                       actionButton("E3", img(src="MfuzzE-3.png", height = 125, width = 125)),
                                       bsModal("e3", "Mfuzz3", "", size = "large", img(src="MfuzzE-3.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="3-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE3')),
                                       actionButton("E4", img(src="MfuzzE-4.png", height = 125, width = 125)),
                                       bsModal("e4", "Mfuzz4", "", size = "large", img(src="MfuzzE-4.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="4-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE4')),
                                       actionButton("E5", img(src="MfuzzE-5.png", height = 125, width = 125)),
                                       bsModal("e5", "Mfuzz5", "", size = "large", img(src="MfuzzE-5.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="5-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE5')),
                                       actionButton("E6", img(src="MfuzzE-6.png", height = 125, width = 125)),
                                       bsModal("e6", "Mfuzz6", "", size = "large", img(src="MfuzzE-6.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="6-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE6')),
                                       actionButton("E7", img(src="MfuzzE-7.png", height = 125, width = 125)),
                                       bsModal("e7", "Mfuzz7", "", size = "large", img(src="MfuzzE-7.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="7-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE7')),
                                       actionButton("E8", img(src="MfuzzE-8.png", height = 125, width = 125)),
                                       bsModal("e8", "Mfuzz8", "", size = "large", img(src="MfuzzE-8.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="8-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableE8'))
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
                              tabPanel("About",
                                       mainPanel(
                                         h4("About this site:"),
                                         p("Blah-Blah-Blah-Blah-Blah-Blah-Blah-Blah-")
                                       )
                              )
                  )),
                mainPanel(
                  h2("Results:"),
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
                              tabPanel("PubMed Hits", 
                                       h4("PubMed Hits"),
                                       DT::dataTableOutput("tableP")
                                      )
                )
                )
                
)
