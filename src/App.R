
# input is fixed to uploading one results file exported from Quantstudio 5, .txt format
# ! will have to check naming of the standards and controls!!!
# output files are in a .zip
# pic of the standard curve
# Excel workbook with all processed info (PMR for all targets, Ct+SD, DNA input, analysis reps)
# Final csv file that is imported for commercial WIDqEC into Swisslab

# load libraries
library(shiny)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)

# source main calculate_PMR function
source("./calculate_pmr_widqec.R")

# Define UI ----
ui <- fluidPage(
  # App title ----
  titlePanel(h3("Calculate PMR and WID-qEC")),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      fileInput("file1", "Choose Results Analysis File", # Input: Select a file ----
                multiple = FALSE,
                accept = c(".txt")),
      tableOutput("resultsQuantstudio"), # for display uploaded data
      tags$hr(), # Horizontal line ----
      downloadButton("download", "Download PMR Result File") #Download button
    ),

    # Main panel for displaying part of output ----
    mainPanel(
      imageOutput("stdCurve") # display standard curve
    ))
)

# Define server logic ----
server <- function(input, output) {
  #Calculate PMR and standard curve upon uploading of files exported from QuantStudio
  myResults <- eventReactive({input$file1}, {
    # avoid displaying error message
    req(input$file1)
    
    #read in inputfiles
    data <- read.table(file = input$file1$datapath,
                        sep = "\t",
                        skip = 47,
                        header = TRUE)
    
    # ensure sample IDs are in correct format
    data$Sample.Name <- as.character(data$Sample.Name)
    
    # calculate PMR
    results <- calculate_pmr(data)
    results
  })
  
  # display uploaded data
  output$resultsQuantstudio <- renderTable({
    req(input$file1)
    data <- read.table(file = input$file1$datapath,
                       sep = "\t",
                       skip = 47,
                       header = TRUE)
    head(data)[,1:6]
  })
  
  # display calculated curve
  output$stdCurve <- renderPlot({
    myResults()[[1]]
  })
  
  # download batch results
  output$download <- downloadHandler(
    
    filename = function() {paste0("Final_results_", Sys.Date(), ".zip")},
    content = function(file){
      
      #locate temporary dir created in this R session & create sub-dir with time within
      temp_dir <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_dir)
      
      # save cal curve to temporary dir
      ggsave(myResults()[[1]],
             file = paste0(temp_dir, "/","calibration_curve.png"),
             width = 5,
             height = 4) 
      
      # save CT plot to temporary dir
      ggsave(myResults()[[10]],
             file = paste0(temp_dir, "/","CTCOL2A1_MeanStdev.png"),
             width = 85,
             height = 80,
             unit = "mm") 
      
      # save batched results to temporary dir
      wb <- createWorkbook()
      
      addWorksheet(wb, "PMR Values")
      writeDataTable(wb = wb, sheet = 1, x = myResults()[[2]], rowNames=FALSE) #PMRs
      addWorksheet(wb, "Cq mean")
      writeDataTable(wb = wb, sheet = 2, x = myResults()[[3]], rowNames=FALSE) #mean Cq
      addWorksheet(wb, "Cq SD")
      writeDataTable(wb = wb, sheet = 3, x = myResults()[[4]], rowNames=FALSE) #Cq stdev
      addWorksheet(wb, "Input conc")
      writeDataTable(wb=wb, sheet = 4, x = myResults()[[5]], rowNames=FALSE) #log(copy number/5uL)
      addWorksheet(wb, "low DNA input")
      writeDataTable(wb=wb, sheet = 5, x = myResults()[[6]], rowNames=FALSE) #samples for which COL2A1 failed in both reps, should have negative controls
      addWorksheet(wb, "Reprocessing needed")
      writeDataTable(wb = wb, sheet = 6, x = myResults()[[7]], rowNames=FALSE) # samples for which for only one of two reps COL2A1 failed
      addWorksheet(wb, "Reprocessing recommended")
      writeDataTable(wb = wb, sheet =7, x = myResults()[[8]], rowNames=FALSE) # samples for which for only one of two reps target amplified
      addWorksheet(wb, "Warning COL2A1 SD high")
      writeDataTable(wb=wb, sheet = 8, x= myResults()[[9]], rowNames=FALSE)# samples for which SD CT COL2A1 > 1.5

      saveWorkbook(wb, paste0(temp_dir, "/","batch_results.xlsx"))
      
      # save final result summary separately to temporary dir
      write.csv(myResults()[[11]], 
                file=paste0(temp_dir, "/", "final_results.csv"), 
                row.names=FALSE, 
                fileEncoding = "UTF-8")
      
      # zip up temporary dir to file
      zip::zip(
        zipfile = file,
        files = dir(temp_dir),
        root = temp_dir
      )
    },
    contentType = "application/zip"
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)