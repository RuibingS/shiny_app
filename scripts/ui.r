## First specify the packages of interest
packages = c("shiny", "ggplot2",
             "scales", "tidyverse","readxl")
    
## Now load or install&loa        d all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


#library(shiny)
#library(ggplot2)
#library(scales)

#setwd("/srv/shiny-server/r_workspace_191216/scripts")
#setwd("/home/rsh19/Downloads/r_shiny/r_workspace/scripts")
source("assayplot.r")
source("concentration.r")
source("data.r")
source("estimator.r")
source("model.r")
source("sigmoid.r")

#input.dataset <- data.list("../db")

# input.dataset <- c("Procymidone" = "xy.procymidone", "Vinclozolin" = "xy.vinclozolin")

help.text <-
  "Lorem ipsum dolor sit amet, consetetur sadipscing elitr,
sed diam nonumy eirmod tempor invidunt ut labore et dolore
magna aliquyam erat, sed diam voluptua. At vero eos et
accusam et justo duo dolores et ea rebum. Stet clita kasd
gubergren, no sea takimata sanctus est Lorem ipsum dolor
sit amet. Lorem ipsum dolor sit amet, consetetur sadipscing
elitr, sed diam nonumy eirmod tempor invidunt ut labore et
dolore magna aliquyam erat, sed diam voluptua. At vero eos
et accusam et justo duo dolores et ea rebum. Stet clita
kasd gubergren, no sea takimata sanctus est Lorem ipsum
dolor sit amet."


input.models <-
  c(
    "[None]" = "none",
    "[Auto]" = "auto",
    "Probit" = "probit",
    "Logit" = "logit",
    "Weibull" = "weibull",
    "Generalized Logit I" = "glogitI",
    "Generalized Logit II" = "glogitII",
    "Aranda-Ordaz" = "ao",
    "Gompertz" = "gompertz"
  )

shinyUI(fluidPage(theme = "bootstrap.css",
                  tags$header(tags$img(src="logo.jpg", style="margin:0.10% 0% 0.10% 0%;height:50px;"), style="background-color:#003366;"), 
                  
                  navbarPage(
                    "Menu",
                    
                    tabPanel(
                      "Assay",
                      #titlePanel("Inhibition assay operation"),

                      
                      sidebarLayout(
                        sidebarPanel(
                          
                          # Input: Select a file ----
                          fileInput("file1", "Choose File"
                                    ),
                          tags$hr(),
                          
                     #     selectInput(
                    #        inputId = "dataset",
                    #        label = "Data set:",
                    #       data.list("../database"), 
                    #        #data.list("../db"),
                    #        multiple = FALSE,
                    #        selectize = TRUE
                    #      ),
                          selectInput(
                            inputId = "compound",
                            label = "Compound:",
                            c(),
                            multiple = FALSE,
                            selectize = TRUE
                          ),
                          selectInput(
                            inputId = "model",
                            label = "Model:",
                            input.models,
                            multiple = FALSE,
                            selectize = TRUE
                          ),
                          tags$hr(),
                          #br(),
                          radioButtons(
                            inputId = "method",
                            label = "Fitting method:",
                            c(
                              "Absolute" = "abs",
                              "Least squares" = "lsq",
                              "Huber loss (robust)" = "huber",
                              "Tukey's biweight (robust)" = "tukey"
                            )
                          ),
                          tags$hr(),
                          checkboxInput(
                            inputId = "sub_blank",
                            label = strong("Substract blank compound (if available)"),
                            value = FALSE
                          ),
                          #tags$hr(),
                          div(style="display:inline-block;margin-left: 70%;padding-bottom: 10px;",
                          actionButton("Button", "Run")),
                          
                          #tags$hr(),
                          textInput(inputId = "x_axis",
                                    label = "x axis",
                                    value = paste0("conc [","\u03BC", "M]")),
                          textInput(inputId = "y_axis",
                                    label = "y axis:",
                                    value = "growth OD600"),
                          #actionButton("Button", "Run"),
                          #br(),
                          tags$hr(),
                          fluidRow(column(6,
                                          checkboxInput(
                                            inputId = "show_data",
                                            label = strong("Show data"),
                                            value = TRUE
                                          ),
                                          checkboxInput(
                                            inputId = "show_curve",
                                            label = strong("Show curve"),
                                            value = TRUE
                                          ),
                                          tags$hr(),
                                          checkboxInput(
                                            inputId = "show_sd",
                                            label = strong("Show standard deviation"),
                                            value = FALSE
                                          ),
                                          checkboxInput(
                                            inputId = "show_band",
                                            label = strong("Show confidence band (slow)"),
                                            value = FALSE
                                          )
                                          # checkboxInput(
                                          #   inputId = "show_X",
                                          #   label = strong("Show TODO"),
                                          #   value = FALSE
                                          # )
                          ),
                          column(6,
                                 checkboxInput(
                                   inputId = "show_ic50",
                                   label = strong("Show IC50"),
                                   value = FALSE
                                 ),
                                 checkboxInput(
                                   inputId = "show_ic90",
                                   label = strong("Show IC90"),
                                   value = FALSE
                                 ),
                                 checkboxInput(
                                   inputId = "show_mic",
                                   label = strong("Show MIC"),
                                   value = FALSE
                                 ),
                                 checkboxInput(
                                   inputId = "show_nic",
                                   label = strong("Show NIC"),
                                   value = FALSE
                                 ))),
                          tags$hr(),
                          
                          numericInput("fheight", "Height (cm)", min=2, max=15, step=1, value = 10),
                          numericInput("fwidth", "Width (cm)", min=2, max=15, step=1, value = 12),
                          selectInput("fres", "Res", choices=c("100","200","300"), selected = "100"),
                          selectInput("fformat", "File type", choices=c("png","tiff","jpeg","pdf"), selected = "png", multiple = FALSE, selectize = TRUE),
                          downloadButton("download", "Download")
                         
                        ),
                        
                        # Show a tabset that includes a plot and table view
                        mainPanel(tabsetPanel(
                          type = "tabs",
                          tabPanel(
                            "Plot",
                            br(),
                            plotOutput(
                              outputId = "plot"
                              #             ,
                              #             height = "500px",
                              #             width = "750px"
                            )
                          ),
                          tabPanel("Summary", br(), verbatimTextOutput("summary"))
                        ))
                      )
                    ),
                    

                    
                    tabPanel(
                      "Help",
                      titlePanel("Help section"),
                      navlistPanel(
                        "Topic A",
                        tabPanel("Section 1", help.text),
                        tabPanel("Section 2", help.text),
                        "Topic B",
                        tabPanel("Section 3", help.text),
                        tabPanel("Section 4", help.text),
                        tabPanel("Section 5", help.text)
                      )
                    )
                    
                  )))
