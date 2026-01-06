library(shiny)

# Define UI
ui <- fluidPage(
  titlePanel("Data Ingestion App"),
  sidebarLayout(
    sidebarPanel(
      textInput("expName", "Experiment Name"),
      fileInput("dataFile", "Upload Data CSV", accept = ".csv"),
      fileInput("ddFile", "Upload Data Dictionary CSV", accept = ".csv"),
      checkboxInput("appendFLMD", "Append FLMD", value = FALSE),
      actionButton("process", "Process Data")
    ),
    mainPanel(
      verbatimTextOutput("status")
    )
  )
)

# Define Server
server <- function(input, output) {
  observeEvent(input$process, {
    req(input$dataFile, input$ddFile, input$expName)
    
    dataPath <- input$dataFile$datapath
    ddPath <- input$ddFile$datapath
    appendFLMD <- input$appendFLMD
    
    output$status <- renderText({
      tryCatch({
        inputDat.fx(input$expName, dataPath, ddPath, appendFLMD)
        "Data processing completed successfully."
      }, error = function(e) {
        paste("Error:", e$message)
      })
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
