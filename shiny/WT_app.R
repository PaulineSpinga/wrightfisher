#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(dplyr)
library(readr)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel(title=h1("Wright Fisher", align="center")),
    
    
    # Sidebar with a slider input for parameters 
    sidebarLayout(
        sidebarPanel(
            numericInput("pop", label = ("Population size"), value = 10, min=0),
            numericInput("allele", label = ("Alleles number"), value = 10,min=0),
            sliderInput("nbG","Generations number:", pre="Value = ",min = 1,max = 1000,value = 10, animate = animationOptions(interval = 50, loop=TRUE)),
            radioButtons("type", "Select the type of graph", choices=list("CoeffH","A","Both"), selected="Anumber"), #, selected="Both"
            actionButton("run", "Display plot")
        ),
        
        
        # Show plots
        mainPanel(
            fluidRow(
                plotOutput("plot"),
                textOutput(outputId = "descPlot"),
                plotOutput("plot2"),
                textOutput(outputId = "descPlot2"),
            ),
        )
    )
)

# Define server logic required to draw plots
server <- function(input, output) {
    
    
    gen <- reactive({
        popSimu(input$pop,input$allele)
    }
    )
    
    
    observe({
        observeEvent(input$run, {
            if (input$type == "A") {
                output$plot <- renderPlot({
                    pA <- plotA(gen())
                },height = 400, width = 700)
                output$descPlot <- renderText({
                    paste("This plot shows alleles number over time for an initial population of ", input$pop, "individuals and ", input$allele, "studied alleles.")
                })
                
                output$plot2 <- NULL
                output$descPlot2 <- NULL
            }
            
            else if (input$type == "CoeffH") {
                output$plot <- renderPlot({
                    gen()
                    pH <- plotH(gen())
                },height = 400, width = 700)
                output$descPlot <- renderText({
                    paste("This plot shows the Heterozygosity Ratio over time for an initial population of ", input$pop, "individuals and ", input$allele, "studied alleles.")
                })
                
                output$plot2 <- NULL
                output$descPlot2 <- NULL
            }
            
            else if (input$type == "Both") {
                output$plot <- renderPlot({
                    pA <- plotA(gen())
                },height = 400, width = 700)
                output$descPlot <- renderText({
                    paste("This plot shows alleles number over time for an initial population of ", input$pop, "individuals and ", input$allele, "studied alleles.")
                })
                
                output$plot2 <- renderPlot({
                    gen()
                    pH <- plotH(gen())
                },height = 400, width = 700)
                output$descPlot2 <- renderText({
                    paste("This plot shows the Heterozygosity Ratio over time for an initial population of ", input$pop, "individuals and ", input$allele, "studied alleles.")
                })
                
                
            }
        })
    })
    
    

    
}

# Run the application 
shinyApp(ui = ui, server = server)
