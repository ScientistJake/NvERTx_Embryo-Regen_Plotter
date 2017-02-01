###NvERTx Shiny Plotter!

library(ggplot2)
library(reshape2)
library(plyr)
library(shinyBS)
library(shiny)

ui <- fluidPage(
  headerPanel('NvERTx Embryo-Regen Plotter'),
  sidebarPanel(
    textInput('search', 'Search for a gene', value = '',width = NULL, placeholder = "e.g. 'TCF'"),
    actionButton("go", "Search!"),
    h4(""),
    textInput('gene1', 'Input NvERTx number', value = "", width = NULL, placeholder =  'e.g. NvERTx.2.133024'),
    textInput('gene2', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
    textInput('gene3', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
    textInput('gene4', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
    textInput('gene5', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
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
    bsModal("modal", "Search Results", "", size = "large", tableOutput('tableA'))
  ),
  mainPanel(
    img(src="Nematostella.png", height = 250, width = 250),
    h4("Explore Mfuzz Clusters:"),
    actionButton("M1", img(src="Mfuzz1.png", height = 125, width = 125)),
    bsModal("m1", "Mfuzz1", "", size = "large", img(src="Mfuzz1.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="1-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM1')),
    actionButton("M2", img(src="Mfuzz2.png", height = 125, width = 125)),
    bsModal("m2", "Mfuzz2", "", size = "large", img(src="Mfuzz2.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="2-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM2')),
    actionButton("M3", img(src="Mfuzz3.png", height = 125, width = 125)),
    bsModal("m3", "Mfuzz3", "", size = "large", img(src="Mfuzz3.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="3-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM3')),
    actionButton("M4", img(src="Mfuzz4.png", height = 125, width = 125)),
    bsModal("m4", "Mfuzz4", "", size = "large", img(src="Mfuzz4.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="4-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM4')),
    actionButton("M5", img(src="Mfuzz5.png", height = 125, width = 125)),
    bsModal("m5", "Mfuzz5", "", size = "large", img(src="Mfuzz5.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="5-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM5')),
    actionButton("M6", img(src="Mfuzz6.png", height = 125, width = 125)),
    bsModal("m6", "Mfuzz6", "", size = "large", img(src="Mfuzz6.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="6-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM6')),
    actionButton("M7", img(src="Mfuzz7.png", height = 125, width = 125)),
    bsModal("m7", "Mfuzz7", "", size = "large", img(src="Mfuzz7.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="7-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM7')),
    actionButton("M8", img(src="Mfuzz8.png", height = 125, width = 125)),
    bsModal("m8", "Mfuzz8", "", size = "large", img(src="Mfuzz8.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="8-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM8')),
    actionButton("M9", img(src="Mfuzz9.png", height = 125, width = 125)),
    bsModal("m9", "Mfuzz9", "", size = "large", img(src="Mfuzz9.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="9-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM9')),
    actionButton("M10", img(src="Mfuzz10.png", height = 125, width = 125)),
    bsModal("m10", "Mfuzz10", "", size = "large", img(src="Mfuzz10.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="10-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM10')),
    actionButton("M11", img(src="Mfuzz11.png", height = 125, width = 125)),
    bsModal("m11", "Mfuzz11", "", size = "large", img(src="Mfuzz11.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="11-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM11')),
    actionButton("M12", img(src="Mfuzz12.png", height = 125, width = 125)),
    bsModal("m12", "Mfuzz12", "", size = "large", img(src="Mfuzz12.png", height = 350, width = 350),h4("Biological process GO enrichment:"),img(src="12-BP-plot-2.png", height = 450, width = 600),tableOutput('tableM12')),
    h2("Regeneration Expression"),
    plotOutput("plot1"),
    h2("Embryonic Expression"),
    plotOutput("plot2"),
    h4("Regeneration average counts (hours post amputation)"),
    tableOutput('table'),
    h4("Warner et al. average counts (hours post fertilization)"),
    tableOutput('table3'),
    h4("Fischer et al. counts (hours post fertilization)"),
    tableOutput('table4'),
    h4("Helm et al. counts (hours post fertilization)"),
    tableOutput('table5'),
    h4("Annotation"),
    tableOutput('table2')
    )
)

