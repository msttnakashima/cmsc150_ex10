#specified behavior of the app 
source("NakashimaEx08.r") #for quadratic spline interpolation
source("NakashimaEx09.r") #for simplex method 

#DIVOC initial values
simplex_divoc_input = matrix(0, nrow = 4, ncol = 6, byrow = TRUE)

rownames = c("Demands", "Denver", "Phoenix", "Dallas")
colnames = c("Supply", "SAC", "SL", "ALB", "CHI", "NYC")

colnames(simplex_divoc_input) = colnames
rownames(simplex_divoc_input) = rownames

simplex_divoc_input[1,] = c(0, 180, 80, 200, 160, 220)
simplex_divoc_input[2,] = c(310, 10, 8, 6, 5, 4)
simplex_divoc_input[3,] = c(260, 6, 5, 4, 3, 6)
simplex_divoc_input[4,] = c(280, 3, 4, 5, 5, 9)

#DIVOC output
simplex_divoc_output = matrix(0, nrow = 4, ncol = 6, byrow = TRUE) 

#GENERIC SOLVER COLNAMES 

#colnames for maximization (x variables first before slack variables)
colnames_gen_max_input <- function (rows, cols) {
  colnames = c()
  for (i in 1:cols) {
    #((number of cols) - (number of rows)) - 1 = number of x variables excluding slack variables & RHS
    if (i <= ((cols - rows) - 1)) colnames = c(colnames, paste("x", i, sep = ""))  #x variables
    else if (i <= (cols - 2)) {  #slack variables
      number = (i - (cols - rows)) + 1
      colnames = c(colnames, paste("S", number, sep = ""))
    }
    else if (i == (cols - 1)) colnames = c(colnames, "Z") #Z column
    else if (i == cols) colnames = c(colnames, "RHS") #RHS column
  }
  return(colnames)
}

#colnames for minimization (slack variables first before x variables)
colnames_gen_min_input <- function (rows, cols) {
  colnames = c()
  for (i in 1:cols) {
    #((number of cols) - (number of rows)) - 1 = number of slack variables excluding x variables & RHS
    if (i <= ((cols - rows) - 1)) colnames = c(colnames, paste("S", i, sep = ""))  #slack variables
    else if (i <= (cols - 2)) {  #x variables
      number = (i - (cols - rows)) + 1
      colnames = c(colnames, paste("x", number, sep = ""))
    }
    else if (i == (cols - 1)) colnames = c(colnames, "Z") #Z column
    else if (i == cols) colnames = c(colnames, "RHS") #RHS column
  }
  return(colnames)
}

#add slack variables 
slack <- function (matrix) {
  for (i in (ncol(matrix) - nrow(matrix)):(ncol(matrix) - 1)) {
    slack = c()
    for (j in 1:nrow(matrix)) {
      slack = c(slack, 0)
    }
    
    index = (i - (ncol(matrix) - nrow(matrix))) + 1
    slack[index] = 1
    matrix[,i] = slack
  }
  return(matrix)
}


server <- function(input, output, session) {
  #QUADRATIC SPLINE INTERPOLATION
  values <- reactiveValues(answer = 0)  #reactive values
  
  #changes in numericInput (input for number of rows)
  observe({
    if (!is.na(input$qsi_row)) {
      temp = matrix(0, nrow = input$qsi_row, ncol = 2, byrow = TRUE)
      colnames(temp) = c("x", "y")
      values$data <- as.data.frame(temp)
    }
  })
  
  #changes in input values in table for qsi
  observeEvent(input$qsi_input_table$changes$changes,{
    values$data = hot_to_r(input$qsi_input_table)
  })
  
  #render input table for qsi 
  output$qsi_input_table <- renderRHandsontable({
    rhandsontable(values$data)
  })
  
  #output qsi functions per interval 
  output$qsifunctions <- renderTable({
    if (sum(values$data$x) == 0 || sum(values$data$y) == 0) return ("Insufficient Data")
    else if (is.na(input$x_qsi)) return ("No value to estimate")
    independentVector = c(values$data$x)
    dependentVector = c(values$data$y)
    data = list(independentVector, dependentVector)
    x_qsi = input$x_qsi 
    
    qsi = poly.qsi(data, x_qsi)
    if (!is.list(qsi)) return (qsi)
    functions_list = qsi$qsi.fxns
    values$answer = qsi$y
    functions_vector = c()
    intervals = c()
    
    #get intervals and corresponding functions per interval
    for (i in 1:length(functions_list) + 1) {
      temp = ""
      temp = paste(temp, values$data$x[i - 1], sep = "")
      temp = paste(temp, " <= x <= ", sep = "")
      temp = paste(temp, values$data$x[i], sep = "")
      intervals[i - 1] = temp
      current = functions_list[[i - 1]]
      functions_vector[i - 1] = deparse(current, width.cutoff = 500)[2]
    }
    result = data.frame(Interval = intervals, Functions = functions_vector)
    return(result)
  }, bordered = TRUE, align = "c")
  
  #output estimated value 
  output$qs_estimate <- renderValueBox({
    out_of_range = valueBox("Error","Value to be estimated is out of range", icon = icon("exclamation-triangle"), color="red", width = 7)
    no_data = valueBox("Error","Insufficient Data", icon = icon("exclamation-triangle"), color="red", width = 7)
    no_estimate = valueBox("Error", "No value to estimate", icon = icon("exclamation-triangle"), color="red", width = 7)
    
    if (sum(values$data$x) == 0 || sum(values$data$y) == 0) return (no_data)
    else if (is.na(input$x_qsi)) return (no_estimate)
    
    independentVector = c(values$data$x)
    dependentVector = c(values$data$y)
    data = list(independentVector, dependentVector)
    x_qsi = input$x_qsi 
    
    qsi = poly.qsi(data, x_qsi)
    if (!is.list(qsi)) return (out_of_range)
    values$answer = qsi$y
    final = valueBox(values$answer, "Estimated Value", icon = icon("slack-hash"), width = NULL)
    return (final)
  })
  
  #SIMPLEX METHOD 
  
  #DIVOC Shipping Analysis
  simplex_values <- reactiveValues(divoc_input = simplex_divoc_input, final.tableau = NULL, basic.solution = NULL, opt.val = NULL, shipping.num = NULL)
  
  #changes in input values in table for simplex
  observeEvent(input$simplex_divoc_input$changes$changes,{
    simplex_values$divoc_input = hot_to_r(input$simplex_divoc_input)
  })
  
  #table input for divoc shipping analysis 
  output$simplex_divoc_input = renderRHandsontable({
    #row & column names
    plants = c("Demands", "Denver","Phoenix","Dallas")
    warehouses = c("Supply", "SAC", "SL", "ALB", "CHI", "NYC")
    data = matrix(simplex_values$divoc_input, nrow = 4, ncol = 6)
    rownames(data) = plants
    colnames(data) = warehouses
    a = rhandsontable(data)
    return(a)
  })
  
  #get minimum cost 
  output$divoc_cost <- renderValueBox({
    no_sol = valueBox("Error","No Feasible Solution",icon = icon("exclamation-triangle"), color="red",width = 7)
    
    simplex_values$opt.val = -1
    
    ship_cost = c()
    
    #get ship costs
    for (i in 2:4) {
      ship_cost = c(ship_cost, simplex_values$divoc_input[i, 2:6])
    }
    
    simplex_values$divoc_ship_cost = ship_cost
    
    #get supply 
    supply = as.vector(simplex_values$divoc_input[2:4,1])
    simplex_values$divoc_supply = supply
    
    #get demands 
    demand = as.vector(simplex_values$divoc_input[1,2:6])
    simplex_values$divoc_demand = demand
    
    #get total demands and supplies
    total_supply = sum(supply)
    total_demand = sum(demand)
    
    
    if (total_supply < total_demand) simplex_values$opt.val = -1
    
    tableau = matrix(nrow = 16, ncol = 25, byrow = TRUE)
    
    #initial tableau for divoc shipping analysis 
    tableau[,1] = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, (demand[1] * -1))  #demand constraint for SAC (flip from >= to <=)
    tableau[,2] = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, (demand[2] * -1))  #demand constraint for SL (flip from >= to <=)
    tableau[,3] = c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, (demand[3] * -1))  #demand constraint for ALB (flip from >= to <=)
    tableau[,4] = c(0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, (demand[4] * -1))  #demand constraint for CHI (flip from >= to <=)
    tableau[,5] = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, (demand[5] * -1))  #demand constraint for NYC (flip from >= to <=)
    tableau[,6] = c(-1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, supply[1])  #supply constraint for DEN (<=)
    tableau[,7] = c(0, 0, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, supply[2])  #supply constraint for PHO (<=)
    tableau[,8] = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, supply[3])  #supply constraint for DAL (<=)
    tableau[,9] = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,10] = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,11] = c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,12] = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,13] = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,14] = c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,15] = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,16] = c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)  #slack variable
    tableau[,17] = c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)  #slack variable
    tableau[,18] = c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)  #slack variable
    tableau[,19] = c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)  #slack variable
    tableau[,20] = c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0)  #slack variable
    tableau[,21] = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)  #slack variable
    tableau[,22] = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)  #slack variable
    tableau[,23] = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0)  #slack variable
    tableau[,24] = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)  #slack variable
    tableau[,25] = c(ship_cost[1], ship_cost[2], ship_cost[3], ship_cost[4], ship_cost[5], ship_cost[6], ship_cost[7], ship_cost[8], ship_cost[9], ship_cost[10], ship_cost[11], ship_cost[12], ship_cost[13], ship_cost[14], ship_cost[15], 0)
    
    simplex_result = simplex(tableau, FALSE, TRUE)
    
    if (!is.list(simplex_result)) return (no_sol)
    simplex_values$final.tableau = simplex_result$final.tableau
    simplex_values$basic.solution = simplex_result$basic.solution
    simplex_values$opt.val = simplex_result$opt.val
    simplex_values$shipping.num = simplex_result$shipping.num
    
    
    
    
    #if total demand is greater than total supply, there's no feasible solution
    if (simplex_values$opt.val == -1) return (no_sol)
    final = valueBox(simplex_values$opt.val, "Minimum Cost", icon = icon("dollar-sign"), width = NULL)
    return (final)
  })
  
  #get number of items shipped per warehouse
  output$divoc_shipping.num = renderRHandsontable({
    if (simplex_values$opt.val == -1) {
      plants = c("DEN","PHO","DAL")
      warehouses = c("SAC", "SL", "ALB", "CHI", "NYC")
      temp = matrix(0, nrow = nrow(simplex_values$shipping.num), ncol = ncol(simplex_values$shipping.num))
      rownames(temp) = plants
      colnames(temp) = warehouses
      final = rhandsontable(temp, readOnly = TRUE)
    }
    else final = rhandsontable(simplex_values$shipping.num, readOnly = TRUE)
    return(final)
  })
  
  #GENERIC SOLVER (GENERAL MAXIMIZATION AND MINIMIZATION PROBLEMS)
  #DIVOC Shipping Analysis
  gen_values <- reactiveValues(opt.val = NULL, gen_basic.solution = NULL, gen_final.tableau = NULL)
  
  #changes in the input for the number of rows and columns
  observe({
    if (!is.na(input$simplex_col) & !is.na(input$simplex_row)) {
      temp = matrix(0, nrow = input$simplex_row, ncol = (input$simplex_col + input$simplex_row + 1), byrow = TRUE)
      if (input$radio_btn == "max") {
        colnames(temp) = colnames_gen_max_input(nrow(temp), ncol(temp)) #if selected radio button is max, call function to generate colnames for max problem
        temp = slack(temp)  #add slack variables (position 1's)
      }
      else {
        colnames(temp) = colnames_gen_min_input(nrow(temp), ncol(temp)) #if selected radio button is min, call function to generate colnames for min problem
        temp = slack(temp) #add slack variables (position 1's)
      }
      gen_values$data <- as.data.frame(temp)
    }
  })
  
  #changes in input values in table for generic solver
  observe({
    if (!is.null(input$simplex_gen_input)) gen_values$data <- hot_to_r(input$simplex_gen_input)
  })
  
  #render input table for generic solver  
  output$simplex_gen_input <- renderRHandsontable({
    rhandsontable(gen_values$data)
  })
  
  #get minimum cost 
  output$opt_val <- renderValueBox({
    gen_values$opt.val = -1
    no_sol = valueBox("Error","No Feasible Solution",icon = icon("exclamation-triangle"), color="red",width = 7)
    
    temp = as.matrix(gen_values$data)
    if (input$radio_btn == "max") isMax = TRUE 
    else isMax = FALSE 
    result = simplex(temp, isMax, FALSE)
    if (!is.list(result)) return(no_sol)
    gen_values$opt.val = result$opt.val
    gen_values$gen_basic.solution = result$basic.solution
    gen_values$gen_final.tableau = result$final.tableau
    
    if (input$radio_btn == "max") final = valueBox(gen_values$opt.val, "Maximum Value", icon = icon("dollar-sign"), width = NULL)
    else final = valueBox(gen_values$opt.val, "Minimum Value", icon = icon("dollar-sign"), width = NULL)
    return (final)
  })
  
  #get basic solution
  output$divoc_basic.solution = renderTable({
    if (input$tabSwitch == "DIVOC Shipping Analysis") {
      if (simplex_values$opt.val == -1) return("No feasible solution")
      else return(simplex_values$basic.solution)
    }
    else if (input$tabSwitch == "Generic Solver") {
      if (gen_values$opt.val == -1) return("No feasible solution")
      else return(gen_values$gen_basic.solution)
    }
  })
  
  #get final tableau
  output$divoc_final.tableau = renderTable({
    if (input$tabSwitch == "DIVOC Shipping Analysis") {
      if (simplex_values$opt.val == -1) return("No feasible solution")
      else return(simplex_values$final.tableau)
    }
    else if (input$tabSwitch == "Generic Solver") {
      if (gen_values$opt.val == -1) return("No feasible solution")
      else return(gen_values$gen_final.tableau)
    }
  })
}