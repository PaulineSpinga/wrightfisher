library(shiny)
library(shinythemes)

shinyUI( fluidPage(
  #theme = bs_theme(version = 4, bootswatch = "journal"),
  #shinythemes::themeSelector(),

  titlePanel(title=h1("Wright Fisher", align="center")),

  sidebarLayout(
    sidebarPanel(
      numericInput("pop", label = ("Population size : (nb_A + nb_a)/2, Tab A-B-C-D"), value = 10, min=0),
      numericInput("allele", label = ("Number of alleles A (nb_A), Tab A-C"), value = 10,min=0),
      sliderInput("nbG","Generations number, Tab A-C", pre="Value = ",min = 1,max = 1000,value = 10, animate = animationOptions(interval = 600, loop=TRUE)),
      sliderInput("bin","Binwidth, Tab A-C",min = 1,max =30,value = 3),
      radioButtons("type", "Select the type of graph, Tab A-C", choices=list("CoeffH","A","Both"), selected="Both"), #, selected="Both"
      actionButton("run", "Display plot"),

      sliderInput("nbS","Simu number, Tab B-D", pre="Value = ",min = 1,max = 1000,value = 100, animate = animationOptions(interval = 600, loop=TRUE)),
      numericInput("step", label = ("Step, Tab B-D"), value = 2, min=1, width=100),
      numericInput("alpha", label = ("Add a selection effect, Tab C-D"), value = 1, min=0, step=0.01)
    ),

    mainPanel(
      tabsetPanel(type="tab",
                  tabPanel("A-GenPop",

                           plotOutput("plot"),
                           textOutput(outputId = "descPlot"),
                           plotOutput("plot2"),
                           textOutput(outputId = "descPlot2"),
                           plotOutput(outputId = "plot3")
                  ),
                  tabPanel("B-Multi GenPop",
                           plotOutput("plotESimuSimu"),
                           plotOutput("plotVSimuSimu")
                  ),
                  tabPanel("C-Select GenPop",
                           plotOutput("plotS"),
                           textOutput(outputId = "descPlotS"),
                           plotOutput("plot2S"),
                           textOutput(outputId = "descPlot2S"),
                           plotOutput(outputId = "plot3S")
                  ),
                  tabPanel("D-Select Multi GenPop",
                           plotOutput("plotESimuSimuS"),
                           plotOutput("plotVSimuSimuS")
                  ),
                  tabPanel("About",
                           uiOutput(outputId = "description"),

                           uiOutput(outputId = "text1", style="color:red"),
                           uiOutput(outputId = "genpop"),
                           uiOutput(outputId = "selectgenpop"),
                           uiOutput(outputId = "popsize"),
                           uiOutput(outputId = "allelesA"),
                           uiOutput(outputId = "generationNumber"),
                           uiOutput(outputId = "binwidth"),

                           uiOutput(outputId = "text2", style="color:red"),
                           uiOutput(outputId = "multigenpop"),
                           uiOutput(outputId = "selectmultigenpop"),
                           uiOutput(outputId = "simuNumber"),
                           uiOutput(outputId = "Step"),

                           uiOutput(outputId = "common", style="color:red"),
                           uiOutput(outputId = "type"),

                           uiOutput("link")
                  )

      )
    )
  )
))
