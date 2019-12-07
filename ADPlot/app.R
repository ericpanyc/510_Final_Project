#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/

library(shiny)
library(limma)
library(Glimma)
library(edgeR)
library(ggfortify)
library(RColorBrewer)



ui <- fluidPage(
    fluidRow(
        column(3, wellPanel(
            selectInput("target", label = h3("Select box"), 
                        choices = list("default" = "default", "AD_condtion" = "AD_condtion", "sex" = "sex"), 
                        selected = "eefault"),
        )),
        column(6,
               plotOutput("plot1", width = 600, height = 300),
        )
    )

   
)
# Define server logic required to draw a histogram
server <- function(input, output) {
    
        output$plot1 <- renderPlot({
            
            if (input$target == "sex") {
                info <- read.csv("./data/foi.txt", header = TRUE, sep = ",", check.names = FALSE)
                table <- read.table("./data/lcpm_table.txt", header = TRUE, sep = ",", check.names = FALSE)
                col <- which(colnames(info) == input$target)
                fact <- as.factor(info[,col])
                full_table <- table
                full_table$sex <- fact
                autoplot(prcomp(table), data=full_table, colour = "sex")
            }
            if (input$target == "AD_condtion") {
                info <- read.csv("./data/foi.txt", header = TRUE, sep = ",", check.names = FALSE)
                table <- read.table("./data/lcpm_table.txt", header = TRUE, sep = ",", check.names = FALSE)
                col <- which(colnames(info) == input$target)
                fact <- as.factor(info[,col])
                full_table <- table
                full_table$AD <- fact
                autoplot(prcomp(table), data=full_table, colour = "AD")
            }
            
            
            
            
        })
    
    
   
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
