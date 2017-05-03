###NvERTx Shiny Plotter!

#It's good practice to have as much code outside the reactive expression as possible.
#I put as much in this top part as possible
library(ggplot2)
library(reshape2)
library(plyr)
library(shinyBS)
library(shinythemes)
library(DT)

#here we start with the shiny
library(shiny)

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
                                  }
                                  
                                  "))
                  ),
                headerPanel('NvERTx Embryo-Regen Plotter'),
                sidebarPanel(
                  tags$head(tags$script(HTML(jscode))),
                  h5("Search for a gene"),
                  tagAppendAttributes(
                    textInput('search', label=NULL, value = '',width = NULL, placeholder = "e.g. 'TCF'"),
                    `data-proxy-click` = "go"
                  ),
                  actionButton("go", "Search!"),
                  h5(a("or click here", href="http://134.59.51.195:4567/")," to blast against NvERTx.2"),
                  h4(""),
                  tags$head(tags$script(HTML(jscode))),
                  tagAppendAttributes(
                    textInput('gene1', 'Input NvERTx number', value = "", width = NULL, placeholder =  'e.g. NvERTx.2.133024'),
                    `data-proxy-click` = "do"
                  ),
                  tags$head(tags$script(HTML(jscode))),
                  tagAppendAttributes(
                    textInput('gene2', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
                    `data-proxy-click` = "do"
                  ),
                  tags$head(tags$script(HTML(jscode))),
                  tagAppendAttributes(
                    textInput('gene3', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
                    `data-proxy-click` = "do"
                  ),
                  tags$head(tags$script(HTML(jscode))),
                  tagAppendAttributes(
                    textInput('gene4', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
                    `data-proxy-click` = "do"
                  ),
                  tags$head(tags$script(HTML(jscode))),
                  tagAppendAttributes(
                    textInput('gene5', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
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
                  )
                ),
                
                mainPanel(
                  bsModal("modal", "Search Results", "", size = "large", DT::dataTableOutput('tableA'))
                ),
                mainPanel(
                  h4("Explore Mfuzz Clusters:"),
                  actionButton("M1", img(src="Mfuzz1.png", height = 125, width = 125)),
                  bsModal("m1", "Mfuzz1", "", size = "large", img(src="Mfuzz1.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="1-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM1')),
                  actionButton("M2", img(src="Mfuzz2.png", height = 125, width = 125)),
                  bsModal("m2", "Mfuzz2", "", size = "large", img(src="Mfuzz2.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="2-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM2')),
                  actionButton("M3", img(src="Mfuzz3.png", height = 125, width = 125)),
                  bsModal("m3", "Mfuzz3", "", size = "large", img(src="Mfuzz3.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="3-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM3')),
                  actionButton("M4", img(src="Mfuzz4.png", height = 125, width = 125)),
                  bsModal("m4", "Mfuzz4", "", size = "large", img(src="Mfuzz4.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="4-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM4')),
                  actionButton("M5", img(src="Mfuzz5.png", height = 125, width = 125)),
                  bsModal("m5", "Mfuzz5", "", size = "large", img(src="Mfuzz5.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="5-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM5')),
                  actionButton("M6", img(src="Mfuzz6.png", height = 125, width = 125)),
                  bsModal("m6", "Mfuzz6", "", size = "large", img(src="Mfuzz6.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="6-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM6')),
                  actionButton("M7", img(src="Mfuzz7.png", height = 125, width = 125)),
                  bsModal("m7", "Mfuzz7", "", size = "large", img(src="Mfuzz7.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="7-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM7')),
                  actionButton("M8", img(src="Mfuzz8.png", height = 125, width = 125)),
                  bsModal("m8", "Mfuzz8", "", size = "large", img(src="Mfuzz8.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="8-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM8')),
                  actionButton("M9", img(src="Mfuzz9.png", height = 125, width = 125)),
                  bsModal("m9", "Mfuzz9", "", size = "large", img(src="Mfuzz9.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="9-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM9')),
                  actionButton("M10", img(src="Mfuzz10.png", height = 125, width = 125)),
                  bsModal("m10", "Mfuzz10", "", size = "large", img(src="Mfuzz10.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="10-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM10')),
                  actionButton("M11", img(src="Mfuzz11.png", height = 125, width = 125)),
                  bsModal("m11", "Mfuzz11", "", size = "large", img(src="Mfuzz11.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="11-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM11')),
                  actionButton("M12", img(src="Mfuzz12.png", height = 125, width = 125)),
                  bsModal("m12", "Mfuzz12", "", size = "large", img(src="Mfuzz12.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="12-BP-plot-2.png", height = 450, width = 600),dataTableOutput('tableM12')),
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
                                       tableOutput('table2')
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
                
))
