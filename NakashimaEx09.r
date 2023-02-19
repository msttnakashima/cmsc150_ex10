#Simplex Method (Minimization & Maximization)
#input is at the bottom of the code

simplex = function(tableau, isMax, problem) {
  tableauRows = nrow(tableau) #number of rows
  tableauCols = ncol(tableau) #number of cols
  
  #GENERATE COLUMN NAMES
  colnames = c()  #vector to store column names 
  
  if (isMax == FALSE) { #minimization: slack variables come before the x variables 
    for (i in 1:tableauCols) {
      #((number of cols) - (number of rows)) - 1 = number of slack variables excluding x variables & RHS
      if (i <= ((tableauCols - tableauRows) - 1)) colnames = c(colnames, paste("S", i, sep = ""))  #slack variables
      else if (i <= (tableauCols - 2)) {  #x variables
        number = (i - (tableauCols - tableauRows)) + 1
        colnames = c(colnames, paste("x", number, sep = ""))
      }
      else if (i == (tableauCols - 1)) colnames = c(colnames, "Z") #Z column
      else if (i == tableauCols) colnames = c(colnames, "RHS") #RHS column
    }
    colnames(tableau) = colnames
  }
  else if (isMax == TRUE) {  #maximization: x variables come before the slack variables
    for (i in 1:tableauCols) {
      #((number of cols) - (number of rows)) - 1 = number of x variables excluding slack variables & RHS
      if (i <= ((tableauCols - tableauRows) - 1)) colnames = c(colnames, paste("x", i, sep = ""))  #x variables
      else if (i <= (tableauCols - 2)) {  #slack variables
        number = (i - (tableauCols - tableauRows)) + 1
        colnames = c(colnames, paste("S", number, sep = ""))
      }
      else if (i == (tableauCols - 1)) colnames = c(colnames, "Z") #Z column
      else if (i == tableauCols) colnames = c(colnames, "RHS") #RHS column
    }
    colnames(tableau) = colnames
  }
  
  #loop until there are no negative values in the bottom row
  while(sum(tableau[tableauRows, -tableauCols] < 0) > 0) {  #tableau[tableauRows, -tableauCols] = bottom row excluding RHS column
    #find the smallest value in the bottom row excluding RHS column
    pivot_column = which.min(tableau[tableauRows, -tableauCols]) 
    smallest = tableau[length(tableau[, 1]), pivot_column] 
    
    column = tableau[-tableauCols,] 
    column = column[, pivot_column]  #store pivot column (without bottom row)
    
    #get test ratio
    testratio = c() #vector to store all test ratios
    for (i in 1:tableauRows - 1) {  #get test ratios 
      testratio[i] = tableau[i, tableauCols] / tableau[i, pivot_column]
    }
    testratio[testratio <= 0] = Inf     #convert negative test ratios to Inf (to not be included in which.min)
    testratio[is.na(testratio)] = Inf   #convert NaN (0/0) values to Inf (to not be included in which.min)
    
    if (all(testratio == testratio[1])) return ("No feasible solution.")  #if all test ratios are Inf, there is no feasible solution
    
    pivot_row = which.min(testratio)   #get smallest positive test ratio
    
    pivot_element = tableau[pivot_row, pivot_column]
    
    #GAUSS-JORDAN ELIMINATION
    
    #find a[i,]: a[i,] / a[i,i]
    tableau[pivot_row,] <- tableau[pivot_row,] / pivot_element
    for(i in 1:tableauRows) {
      if(i == pivot_row) next
      
      #find NORMALIZED ROW: a[j,i] * a[i,]
      normalized_row = tableau[pivot_row,] * tableau[i, pivot_column]
      
      #find a[j,]: a[j,] - NORMALIZED ROW 
      tableau[i,] = tableau[i,] - normalized_row
    }
  }
  
  #BASIC SOLUTION
  solution = c()
  if (isMax == FALSE) { #minimization (refer to bottom row)
    solution = tableau[tableauRows, -(tableauCols - 1)]
  }
  
  else if (isMax == TRUE) {   #maximization (refer to RHS Column)
    #if there is only one 1 in a column, add the corresponding value in the RHS column to "solution" vector
    for(i in 1:(tableauCols - 1)) {
      if (length(unique(tableau[, i])) == 2 & sum(unique(tableau[, i])) == 1) {  #if there is only one 1 in the column
        for (j in 1:nrow(tableau)) {  #find row index of the one 1 in the column
          if (tableau[j, i] == 1) {
            solution = c(solution, tableau[j, tableauCols]) 
            break
          }
        }
      }
      else solution = c(solution, 0) #there is other values in the column besides one 1
    }
  }
  
  basic_sol = matrix(solution, nrow = 1, ncol = length(solution))
  opt.val = basic_sol[1, ncol(basic_sol)] #get optimum value 
  colnames(basic_sol) = colnames[1:(length(colnames) - 1)]
  
  #SHIPPING.NUM (check "problem" to see if "shipping.num" will be returned)
  if (problem == TRUE) {  #if problem is true, return "shipping.num"
    shipping.num = matrix(solution[(tableauCols - tableauRows):(tableauCols - 2)], nrow = 3, ncol = 5, byrow = TRUE)
    rownames(shipping.num) = c("DEN", "PHO", "DAL") #rows represent the plants
    colnames(shipping.num) = c("SAC", "SL", "ALB", "CHI", "NYC")  #cols are the warehouses
    values = list(final.tableau = tableau, basic.solution = basic_sol, opt.val = opt.val, shipping.num = shipping.num)
  }
  else values = list(final.tableau = tableau, basic.solution = basic_sol, opt.val = opt.val) #does not return "shipping.num"
  
  return(values)
}

#USER INPUT 
#NOTE: 
  #For minimization problems, the input (initial tableau) should already be transposed with the slack variables (see sample input below)

# #SAMPLE INPUT: Minimization Problem (initial tableau is already transposed)
# #Constraints:
#    #x1 + 2 * x2 = 4                         1  2  4         1  7  14        1  7   1  0  0  14
#    #7 * x1 + 6 * x2 = 20              =>    7  6 20    =>   2  6  20   =>   2  6   0  1  0  20
#    #MINIMIZE: Z = 14 * x1 + 20 * x2        14 20  1         4  20  1        4  20  0  0  1  0
#                                        #not transposed    #transposed     #add slack variables (final initial tableau / input)
# 
# data = c(1, 7, 1, 0, 0, 14,
#           2, 6, 0, 1, 0, 20,
#           -4, -20, 0, 0, 1, 0
# )
# 
# tableau = matrix(data, nrow = 3, ncol = 6, byrow = TRUE)


isMax = FALSE      #true - max, false - min
problem = TRUE     #true - return shipping.num, false - does NOT return shipping.num


##DIVOC Shipping Analysis 
supply = c(310, 260, 280)
demand = c(180, 80, 200, 160, 220)
ship_cost = c(10, 8, 6, 5, 4, 
              6, 5, 4, 3, 6, 
              3, 4, 5, 5, 9
              )

tableau = matrix(0, ncol = 25, nrow = 16)

#in the initial tableau, negate constraints with ">=" to flip the inequality to "<="
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

result = simplex(tableau, isMax, problem)
result