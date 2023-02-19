library(shiny)
library(shinydashboard)
library(rhandsontable)

ui <- dashboardPage(
  dashboardHeader(title = "CMSC 150"),
  dashboardSidebar(
    tags$style(HTML('.logo, .navbar{background-color: #A1B5D8 !important;}')),
    sidebarMenu(
      menuItem("Quadratic Spline Interpolation", tabName = "qsi", icon = icon("chart-line")),
      menuItem("Simplex Method", tabName = "simplex", icon = icon("sort-numeric-up-alt")),
      menuItem("About", tabName = "about", selected = TRUE)
    )
  ),
  dashboardBody(
    #QUADRATIC SPLINE INTERPOLATION PAGE 
    tabItems(
      tabItem(tabName = "qsi",
        #change colors of boxes
        tags$style(HTML(".box.box-solid.box-primary>.box-header {color:#FFFCF7; background:#738290}
                         .box.box-solid.box-primary {border-bottom-color:#738290;
                                                     border-left-color:#738290;
                                                     border-right-color:#738290;
                                                     border-top-color:#738290;}
        ")),
        
        fluidRow(
          column(width = 2, 
            box(
              width = NULL, 
              status = "primary",
              solidHeader = TRUE, 
              numericInput("qsi_row", "Number of Rows:", min = 1, value = 5),    #input number of rows 
            ),
            
            box(
              width = NULL,
              status = "primary",
              solidHeader = TRUE,
              numericInput("x_qsi", "Value to be estimated:", min = 1, value = 5)    #input value to be approximated
            )
          ),
          
          column(width = 10,
            valueBoxOutput("qs_estimate", width = 12),
            box(
              style = "max-height: 200px; overflow-y: scroll; overflow-x: hidden",
              width = 5,
              title = "Data",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              rHandsontableOutput("qsi_input_table"),     #qsi input table
            ),
            
            box(
              style = "max-height: 200px; overflow-y: scroll;",
              width = 7,
              title = "Polynomial Regression Functions",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              tableOutput("qsifunctions")
            )
          )
        )
      ),
      
      #SIMPLEX METHOD PAGE 
      tabItem(tabName = "simplex",
        fluidRow(
          column(width = 7, 
            tabBox(
              #title = "Input", 
              id = "tabset1",
              width = 12,
              
              tabsetPanel(id = "tabSwitch", 
                #DIVOC SHIPPING ANALYSIS INPUT TAB
                tabPanel("DIVOC Shipping Analysis", 
                  style = "height: 250px; overflow-x: hidden; overflow-y: hidden;",
                  br(),
                  
                  #table input
                  rHandsontableOutput("simplex_divoc_input")
                ),
                
                #MAXIMIZATION & MINIMIZATION INPUT TAB
                tabPanel("Generic Solver",
                  width = 12,
                  id = "generic",
                  style = "height: 250px; overflow-x:hidden; overflow-y: scroll;",
                  fluidRow(
                    column(width = 12,
                      strong("Input initial tableau: "),
                      p(em("Note:"), " For minimization problems, the matrix from the equations should already be transposed."),
                      #br(),
                      
                      rHandsontableOutput("simplex_gen_input"),
                      #style = "height: 250px; overflow-x: scroll;",
                    )
                  )
                )
              )
            )
          ),

          column(width = 5,
            conditionalPanel(
              condition = "input.tabSwitch == 'DIVOC Shipping Analysis'",
                valueBoxOutput("divoc_cost", width = NULL),
                box(
                  width = "NULL",
                  title = "Number of Items Shipped per Warehouse",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  rHandsontableOutput("divoc_shipping.num")
                )
            ),
            
            conditionalPanel(
              condition = "input.tabSwitch == 'Generic Solver'",
              #radio buttons to choose whether to perform minimization or maximization
              valueBoxOutput("opt_val", width = NULL),
              box(
                width = NULL, 
                status = "primary",
                solidHeader = TRUE, 
                
                #option to choose whether to perform maximization or minimization
                radioButtons("radio_btn", "Choose one:",
                             choiceNames = list("Minimization", "Maximization"),
                             choiceValues = list("min", "max"),
                             inline = TRUE
                ),
                
                fluidRow(
                  column(5, 
                    box(
                      width = NULL, 
                      status = "primary",
                      solidHeader = TRUE,
                      numericInput("simplex_row", "No. of Rows:", min = 1, value = 3),    #input number of rows 
                    )
                  ),
                  
                  column(7,
                    box(
                      width = NULL,
                      status = "primary",
                      solidHeader = TRUE,
                      numericInput("simplex_col", "No. of Decision Variables:", min = 2, value = 2),    #input number of columns 
                    )
                  )
                )  
              )
            )
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Basic Solution",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            column(width = 12, 
              tableOutput("divoc_basic.solution"),
              style = "overflow-x: scroll;"
            )
          )
        ),

        fluidRow(
          box(
            width = 12,
            title = "Final Tableau",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            column(width = 12, 
              tableOutput("divoc_final.tableau"),
              style = "overflow-x: scroll"
            )
          )
        )
      ),
      
      #ABOUT PAGE 
      tabItem(tabName = "about",
        tabPanel("About",
          #titlePanel("About"),
          sidebarLayout(
            sidebarPanel(
              h2("Student Info"),
              p(strong("Name: "), "Marie Sophia Therese T. Nakashima"),
              p(strong("Student Number: "), "2020-00822"),
              p(strong("Section: "), "CMSC 150 B-3L"),
              br(),
              br(),
              img(src = "https://pbs.twimg.com/media/FAdpbT5VIAMqQB_?format=jpg&name=small", 
                height = 200, width = "100%",
                style = "object-fit: cover;"
              ),
              p("Credits: ", a("@Togeplzcomeback", href = "https://twitter.com/Togeplzcomeback")),
            ), 
            mainPanel(
              h1("Exercise 10"),
              p("A web application with a user interface for both", 
                strong("Quadratic Spline Interpolation"), 
                "and", 
                strong("Simplex Method"), 
                "written in R with the Rshiny library for GUI."),
              br(), 
              br(),
              h2("Features"), 
              p("- Find polynomials on each interval given a data set using Quadratic Spline Interpolation."),
              p("- Generic solver for maximization problems using Simplex Method."),
              p("- Generic solver for minimization problems using Simplex Method."), 
              p("- Minimize total shipping cost of", 
                em("Dedmond Integrated Valley Operations Company (DIVOC)"),
              "using Simplex Method"),
            )
          )
        )
      )
    )
  )
)