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