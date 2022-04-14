library(parallel)
library(ggplot2)
library(shiny)
library(shinythemes)
library(dplyr)
library(readr)
library(bslib)

library(ggpmisc)
library(expint)
library(wrightfisher)

shinyServer(function(input, output) {

  output$description <- renderText(HTML(paste(strong("In this Shiny app, we study the allele distribution over time of a finite population of size N containing only two types of alleles (A and a) with a dynamics described by the Wright-Fisher model. Our main goal is the study of a quantity called the fixation time."),"<br/>")))

  output$text1 <- renderText(HTML(paste("<br/>",strong("The following parameters allows to study the fixation time for one inital probability"))))
  output$genpop <- renderText(HTML(paste(strong("Tab GenPop :"), "Without Selection")))
  output$selectgenpop <- renderText(HTML(paste(strong("Tab Select GenPop :"), "With Selection")))
  output$popsize <- renderText(HTML(paste(strong("Population size :"), "Number of individuals in the initial population.")))
  output$allelesA <- renderText(HTML(paste(strong("Number of alleles A :"), "Reflects the number of copies of the allele A in the initial population")))
  output$generationNumber <- renderText(HTML(paste(strong("Generations Number:"), "Number of initial populations generated according to the parameters Population size and Number of alleles A")))
  output$binwidth <- renderText(HTML(paste(strong("Binwidth :"), "Bandwidth of the histogram representing the distribution of the fixation time according to the set of initial populations")))
  output$type <- renderText(HTML(paste(strong("Type of graph :")," CoeffH and A allows to view the coefficient of heterozygosity or the number of alleles A over time for the set of initial populations,respectively")))

  output$text2 <- renderText(HTML(paste("<br/>",strong("The following parameters allows to study the fixation time for every possible inital probabilities "))))
  output$multigenpop <- renderText(HTML(paste(strong("Tab Multi GenPop :"), "Without Selection")))
  output$selectmultigenpop <- renderText(HTML(paste(strong("Tab Select Multi GenPop :"), "With Selection")))
  output$simuNumber <- renderText(HTML(paste(strong("Simu Number:"), "Number of initial populations generated according to the parameters Population size and Step")))
  output$Step <- renderText(HTML(paste(strong("Step:"), "Value used to calculate the each initial probability of allele A according to the size of the population")))

  output$common <- renderText(HTML(paste("<br/>",strong("Shared parameters"))))
  output$type <- renderText(HTML(paste(paste(strong("Add a selection effect :"),"Value that increases the probability of transmission of the allele A"),"<br/>")))

  url <- a("GitHub", href="https://github.com/PaulineSpinga/wrightfisher/blob/main/rapport/L3DLSDVINFO-Spinga-Pauline-rapport.pdf")
  output$link <- renderUI({
    tagList("For more details: ", url)
  })

  out <- reactive({
    mclapply(rep(input$pop,input$nbG),popSimu, nbA = input$allele,mc.cores = 1)
  })

  outs <- reactive({
    mclapply(rep(input$pop,input$nbG),popSimu.select, nbA = input$allele, s = input$alpha/(2*input$pop),mc.cores = 1)
  })

  temps_fix <- reactive({vectFix(out())})
  temps_fixs <- reactive({vectFix(outs())})

  output$plotESimuSimu <- renderPlot(
    { withProgress(message = 'Making plot (Fixation time): ', value = 0,{
      valA <- seq(0,2*input$pop, by=input$step)
      if(((2*input$pop)%% input$step)!= 0){
        valA <- c(valA,2*input$pop)
      }
      n <- length(valA)
      outF <- vector("list", length(valA))
      for (i in 1:length(valA)) {
        outF[[i]] <- mclapply(rep(input$pop,input$nbS),popSimu, nbA = valA[i],mc.cores = 1)
        incProgress(1/n, detail = paste(round(i/n*100,0), " %"))
      }
    })
      affichTP(outF, input$pop, input$step)},
    height = 400, width = 700)

  output$plotVSimuSimu <- renderPlot(
    { withProgress(message = 'Making plot (Variance): ', value = 0, {
      valA <- seq(0,2*input$pop, by=input$step)
      if(((2*input$pop)%% input$step)!= 0){
        valA <- c(valA,2*input$pop)
      }
      n <- length(valA)
      outF <- vector("list", length(valA))
      for (i in 1:length(valA)) {
        outF[[i]] <- mclapply(rep(input$pop,input$nbS),popSimu, nbA = valA[i],mc.cores = 1)
        incProgress(1/n, detail = paste(round(i/n*100,0), " %"))
      }
    })
      affichVP(outF)},
    height = 400, width = 700)



  output$plotESimuSimuS <- renderPlot(
    { withProgress(message = 'Making plot: ', value = 0, {
      valA <- seq(0,2*input$pop, by=input$step)
      if(((2*input$pop)%% input$step)!= 0){
        valA <- c(valA,2*input$pop)
      }
      n <- length(valA)
      outF <- vector("list", length(valA))
      outFS <- vector("list", length(valA))
      for (i in 1:length(valA)) {
        outFS[[i]] <- mclapply(rep(input$pop,input$nbS),popSimu.select, nbA = valA[i],s=input$alpha/(2*input$pop),mc.cores = 1)
        #outF[[i]] <- mclapply(rep(input$pop,input$nbS),popSimu, nbA = valA[i],mc.cores = 1)
        incProgress(1/n, detail = paste(round(i/n*100,0), " %"))
      }
    })
      affichTP.both(outF,outFS,input$pop, input$alpha/(2*input$pop))},
    height = 400, width = 700)

  output$plotVSimuSimuS <- renderPlot(
    { withProgress(message = 'Making plot: ', value = 0, {
      valA <- seq(0,2*input$pop, by=input$step)
      if(((2*input$pop)%% input$step)!= 0){
        valA <- c(valA,2*input$pop)
      }
      n <- length(valA)
      outFS <- vector("list", length(valA))
      for (i in 1:length(valA)) {
        outFS[[i]] <- mclapply(rep(input$pop,input$nbS),popSimu.select, nbA = valA[i],s=input$alpha/(2*input$pop),mc.cores = 1)
        incProgress(1/n, detail = paste(round(i/n*100,0), " %"))
      }
    })
      affichVP(outFS)},
    height = 400, width = 700)



  observe({

    observeEvent(input$run, {
      if (input$type == "A") {
        output$plot <- renderPlot({
          affichDataA(out(),input$nbG)
        },height = 400, width = 700)
        output$descPlot <- renderText({
          paste("This plot shows alleles number over time for an initial population of ", input$pop, "individuals and ", input$allele, " alleles A.")
        })

        output$plot2 <- renderPlot({affichDataT(temps_fix(),input$bin)},height = 400, width = 700)
        output$descPlot2 <- NULL

        output$plot3 <- NULL
      }

      else if (input$type == "CoeffH") {
        output$plot <- renderPlot({
          affichDataH(out(),input$nbG)
        },height = 400, width = 700)
        output$descPlot <- renderText({
          paste("This plot shows the Heterozygosity Ratio over time for an initial population of ", input$pop, "individuals and ", input$allele, " alleles A.")
        })

        output$plot2 <- renderPlot({affichDataT(temps_fix(),input$bin)},height = 400, width = 700)
        output$descPlot2 <- NULL

        output$plot3 <- NULL
      }

      else if (input$type == "Both") {
        output$plot <- renderPlot({
          affichDataA(out(),input$nbG)
        },height = 400, width = 700)
        output$descPlot <- renderText({
          paste("This plot shows alleles number over time for an initial population of ", input$pop, "individuals and ", input$allele, " alleles A.")
        })

        output$plot2 <- renderPlot({
          affichDataH(out(), input$nbG)
        },height = 400, width = 700)
        output$descPlot2 <- renderText({
          paste("This plot shows the Heterozygosity Ratio over time for an initial population of ", input$pop, "individuals and ", input$allele, " alleles A.")
        })

        output$plot3 <- renderPlot({affichDataT(temps_fix(),input$bin)},height = 400, width = 700)
      }
    })
  })


  observe({

    observeEvent(input$run, {
      if (input$type == "A") {
        output$plotS <- renderPlot({
          affichDataA(outs(),input$nbG)
        },height = 400, width = 700)
        output$descPlotS <- renderText({
          paste("This plot shows alleles number over time for an initial population of ", input$pop, "individuals and ", input$allele, "alleles A.")
        })

        output$plot2S <- renderPlot({affichDataT(temps_fix(),input$bin)},height = 400, width = 700)
        output$descPlot2S <- NULL

        output$plot3S <- NULL
      }

      else if (input$type == "CoeffH") {
        output$plotS <- renderPlot({
          affichDataH(outs(),input$nbG)
        },height = 400, width = 700)
        output$descPlotS <- renderText({
          paste("This plot shows the Heterozygosity Ratio over time for an initial population of ", input$pop, "individuals and ", input$allele, "alleles A.")
        })

        output$plot2S <- renderPlot({affichDataT(temps_fix(),input$bin)},height = 400, width = 700)
        output$descPlot2S <- NULL

        output$plot3S <- NULL
      }

      else if (input$type == "Both") {
        output$plotS <- renderPlot({
          affichDataA(outs(),input$nbG)
        },height = 400, width = 700)
        output$descPlotS <- renderText({
          paste("This plot shows alleles number over time for an initial population of ", input$pop, "individuals and ", input$allele, "alleles A.")
        })

        output$plot2S <- renderPlot({
          affichDataH(outs(), input$nbG)
        },height = 400, width = 700)
        output$descPlot2S <- renderText({
          paste("This plot shows the Heterozygosity Ratio over time for an initial population of ", input$pop, "individuals and ", input$allele, "alleles A.")
        })

        output$plot3S <- renderPlot({affichDataT(temps_fix(),input$bin)},height = 400, width = 700)
      }
    })
  })




})
