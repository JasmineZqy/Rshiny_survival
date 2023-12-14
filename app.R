library(shiny)
library(survival)
library(flexsurv)
library(survminer)
library(ggplot2)
library(plotly)
library(DT)
source("code.R")

# Define UI for application that draws a histogram
ui <- navbarPage(title ="CardioGuardian",tabPanel("Prediction",fluidPage(
  tags$style(HTML("
    .survtable {
      font-size: 20px;
    }
  ")),
  tags$style("
    .survtabletitle{
      font-size: 20px;
      font-weight: bold;
    }"
  ),
  
  # Application title
  #headerPanel(fluidRow(h1("CardioGuardian: One-stop heart failure prediction and patient information management system")),
  #          fluidRow(h2("Patient Information & Prediction"))),
  fluidRow(
    column(12, align = "left",
           tags$h1("CardioGuardian: One-stop heart failure prediction and clinical data management system"),
           tags$h3("Patient Information Input & Prediction")
    )
  ),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("Name","Name:"),
      sliderInput("Age",
                  "Age:",
                  min = 1,
                  max = 120,
                  value = 50
      ),
      selectInput("Gender","Gender:",choices = list("Male","Female")),
      selectInput("rbc","Decrease of red blood\ncells or hemoglobin:",choices = list("Yes","No")),
      selectInput("hyper","If the patient\nhas hypertension:",choices = list("Yes","No")),
      numericInput("cpk","Level of the CPK enzyme\nin the blood (mcg/L)",value = 0),
      selectInput("diabetes","If the patient\nhas diabetes:",choices = list("Yes","No")),
      sliderInput("bloodleaving","Percentage of blood leaving\nthe heart at each contraction:",min = 0, max =             100,value = 50),
      numericInput("platelets","Platelets in the\nblood (kiloplatelets/mL):",value = 0),
      numericInput("creatinine","Level of serum creatinine\nin the blood (mg/dL):",value = 0),
      numericInput("sodium","Level of serum sodium\nin the blood (mEq/L):",value = 0),
      actionButton("Submit","Submit to Database", style = "background-color: #09f;color: white;")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      h3(textOutput("graphtitle")),
      plotOutput("survPlot"),
      h3(textOutput("probtitle")),
      numericInput("Month","Predicted Month",value = 1, min = 1, max = 12),
      #tableOutput("survPlot_tbl"),
      fluidRow(
        class = "survtabletitle",
        column(2,""),
        column(4,"Individual Survival Prob.",style = ""),
        column(4,"Population Survival Prob."),
        column(2, "Percentile")
      ),
      fluidRow(
        class = "survtable",
        column(2,textOutput("inmonth")),
        column(4,textOutput("patprob")),
        column(4,textOutput("popprob")),
        column(2,textOutput("perc"))
      ),
      plotly::plotlyOutput("survPlot_hist")
      
      
    )
  )
)),
tabPanel("Database",
         fluidPage(fluidRow(
           column(12, align = "left",
                  tags$h1("CardioGuardian: One-stop heart failure prediction and clinical data management system"),
                  tags$h3("Clinical Data Management System")
           )
         ),
         mainPanel(tableOutput("dataset"),
                   actionButton("Remove","Remove last row", style = "background-color: #09f;color: white;"),
                   actionButton("Sort","Sort by Survival",style = "background-color: #09f;color: white;"),
                   actionButton("Clear","Clear", style = "background-color: white;color: red;"),
                   downloadButton("Download","Download the dataset"),
                   fileInput("file","Upload dataset",accept = c("text/csv", "text/comma-separated-values", "text/plain", ".csv"))
         ))))

# Define server logic required to draw a histogram
server <- function(input, output) {
  new_df <- reactive({
    new <- test_df(input$cpk,input$Age,  input$creatinine, input$sodium)
    return(new)
  })
  time<- reactive({
    time<- input$Month*30
    return(time)
  })
  
  output$graphtitle<-renderText({
    
    paste("Survival curve for",input$Name)
  })
  output$inmonth<-renderText({
    paste("In",input$Month,"month(s)")
  })
  output$patprob<-renderText({
    paste0(round(pred(new_df(),time())[[1]]*100,2),"%")
  })
  output$popprob<-renderText({
    paste0(round(pred(new_df(),time())[[3]]*100,2),"%")
  })
  output$perc<-renderText({
    paste0(round(pred(new_df(),time())[[4]]*100,2),"%")
    
  })
  output$survPlot <- renderPlot({
    Makesurvplot(new_df())
  })
  
  output$probtitle<-renderText({
    
    paste("Survival rates of",input$Name,"in future months")
  })
  
  output$survPlot_tbl <- renderTable({
    Makesurvtable(pred(new_df(),time()))
  })
  output$survPlot_hist<-renderPlotly({
    MakeHist(new_df(),time())
    
  })
  ######create databse#####
 rv <- reactiveValues(data = data.frame(Name = NULL,
                                         Age = numeric(0),
                                         Gender = NULL,
                                         RBC = NULL,
                                         hyper = NULL,
                                         cpk = numeric(0),
                                         diabetes = NULL,
                                         pbl = numeric(0),
                                         platelets = numeric(0),
                                         creatinine = numeric(0),
                                         sodium = numeric(0),
                                         Survival = numeric(0)
                                         
 ))
 observeEvent(input$Submit, {
   rv$data <- rbind(rv$data, data.frame(Name = input$Name,
                                        Age = input$Age,
                                        Gender = input$Gender,
                                        RBC = input$rbc,
                                        hyper = input$hyper,
                                        cpk = input$cpk,
                                        diabetes = input$diabetes,
                                        pbl= input$bloodleaving,
                                        platelets = input$platelets,
                                        creatinine = input$creatinine,
                                        sodium = input$sodium,
                                        Survival = round(pred1m(new_df()),time())))
   
   
 })
 
 
  
  
  output$dataset <- renderTable({
    rv$data
  }, rownames = TRUE)
  
  observeEvent(input$Clear, {
    rv$data <- data.frame(Name = NULL,
                          Age = numeric(0),
                          Gender = NULL,
                          RBC = NULL,
                          hyper = NULL,
                          cpk = numeric(0),
                          diabetes = NULL,
                          pbl = numeric(0),
                          platelets = numeric(0),
                          creatinine = numeric(0),
                          sodium = numeric(0),Survival = numeric(0))
  })
  
  observeEvent(input$Remove, {
    
    if (nrow(rv$data) > 0) {
      rv$data <- rv$data[-nrow(rv$data), ]
    }
  })
  
  
  observeEvent(input$Sort, {
    if (nrow(rv$data) > 0) {
      rv$data <- rv$data[order(rv$data$Survival), ]
    }
  })
  output$Download <- downloadHandler(
    filename = function() {
      paste("dataset-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(rv$data, file, row.names = FALSE) 
    }
  )
  observeEvent(input$file, {
    if (!is.null(input$file)) {
      uploaded_data <- read.csv(input$file$datapath)
      
      rv$data <- rbind(rv$data, uploaded_data)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)