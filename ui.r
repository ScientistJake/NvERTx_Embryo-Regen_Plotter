###NvERTx Shiny Plotter!

#It's good practice to have as much code outside the reactive expression as possible.
#I put as much in this top part as possible
library(ggplot2)
library(reshape2)
library(plyr)

#here we start with the shiny
library(shiny)

#this makes the buttons and input
ui <- fluidPage(
  headerPanel('NvERTx Embryo-Regen Plotter'),
  sidebarPanel(
    checkboxInput('log', 'Check for Log2', value = FALSE, width = NULL),
    textInput('gene1', 'Input NvERTx number', value = "", width = NULL, placeholder =  'e.g. NvERTx.2.133024'),
    textInput('gene2', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
    textInput('gene3', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
    textInput('gene4', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
    textInput('gene5', 'Input another number or leave blank', value = "", width = NULL, placeholder = 'e.g. NvERTx.2.133024'),
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
    img(src="Nematostella.png", height = 300, width = 300),
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
