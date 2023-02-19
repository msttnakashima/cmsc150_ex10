#load packages
library(shiny) 
library(shinydashboard)
library(rhandsontable)

#source ui & server
source("ui.r")
source("server.r")

#run app
shinyApp(ui, server)