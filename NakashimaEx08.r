#Quadratic Spline Interpolation 
#input is at the bottom of the code 

poly.qsi = function(data, x) {
  x_values = data[[1]]
  y_values = data[[2]]
  
  if (length(x_values) != length(y_values)) return(NA)  #if x values and y values do not have the same length, return NA
  else if ((x > x_values[length(x_values)]) || (x < x_values[1])) return ("Value to be approximated is out of range")  #if "x" is not within the range of x values, return NA  
  
  n = length(x_values)
  no_of_intervals = n - 1 #no_of_intervals = n - 1
  no_of_equations = no_of_intervals * 3  #3(n) equations
  equations = matrix(0, nrow = no_of_equations, ncol = no_of_equations + 1) #matrix to store coefficients of equations
  internal_and_endpoints = 2 * (n - 2) + 2 #number of internal knots and end points equations
  
  #FORMATION OF EQUATIONS
  #internal knots & end points 
  col = 1
  row = 1
  for (i in seq(1, internal_and_endpoints, by = 2)) { #from 1 to "internal_and_endpoints", increment by 2
    #internal knot: (a[i] * (x[i - 1])^2) + (b[i] * x[i - 1]) + c[i - 1] = f(x[i - 1])
    #end point: (a[1] * (x[0])^2) + (b[1] * x[0]) + c[1] = f(x[0])
    equations[i, col] = x_values[row]^2
    equations[i, col + 1] = x_values[row]
    equations[i, col + 2] = 1
    equations[i, ncol(equations)] = y_values[row]
    
    #internal knot: (a[i - 1] * (x[i - 1])^2) + (b[i - 1] * x[i - 1]) + c[i - 1] = f(x[i - 1])
    #end point: (a[n] * (x[n])^2) + (b[n] * x[n]) + c[n] = f(x[n])
    equations[i + 1, col] = x_values[row + 1]^2
    equations[i + 1, col + 1] = x_values[row + 1]
    equations[i + 1, col + 2] = 1
    equations[i + 1, ncol(equations)] = y_values[row + 1]
    
    col = col + 3   #columns have a pattern of: a[i], b[i], c[i] thus increment by 3
    row = row + 1 
  }
  
  #1st derivative at the interior knots
  col = 1
  for(i in ((internal_and_endpoints + 1):(no_of_equations - 1))) {  
    #1st derivative at the interior knots: (2 * a[i - 1] * x[i - 1]) + b[i - 1] = (2 * a[i]  * x[i - 1]) + b[i]
    equations[i, col] = 2 * x_values[2 + (i - (internal_and_endpoints + 1))] 
    equations[i, col + 1] = 1 
    equations[i, col + 3] = -2 * x_values[2 + (i - (internal_and_endpoints + 1))]  #multiply by negative since it's at the right of the equal sign
    equations[i, col + 4] = -1  #multiply by negative since it's at the right of the equal sign
    col = col + 3   #columns have a pattern of: a[i], b[i], c[i] thus increment by 3
  }
  
  equations[no_of_equations, 1] = 1   #a1 = 0
  
  #GAUSS JORDAN ELIMINATION
  n <- nrow(equations)  #n is equal to the no. of unknown variables (also equal to the no. of rows)
  
  for (i in 1:n) {
    if (i != n) {
      #find PIVOT ROW: row with max(abs(a[i:n,i]))
      pivot_row <- which.max(abs(equations[i:n,i])) + i - 1
      if (equations[pivot_row, i] == 0) return (NA) #no unique solution exists, STOP
      
      #do PARTIAL PIVOTING: swap(PIVOT ROW, a[i,])
      if (pivot_row != i) equations[c(i, pivot_row),] <- equations[c(pivot_row, i),]
    }
    
    #find a[i,]: a[i,] / a[i,i]
    equations[i,] <- equations[i,] / equations[i,i]
    for (j in 1:n) {
      if (i == j) next
      
      #find NORMALIZED ROW: a[j,i] * a[i,]
      normalized_row <- equations[j,i] * equations[i,]
      
      #find a[j,]: a[j,] - NORMALIZED ROW
      equations[j,] <- equations[j,] - normalized_row
    }
  }
  
  #get coefficients
  coefficients_vector <- c()
  for (i in 1:n) {
    coefficients_vector[i] <- equations[i, ncol(equations)]
  }
  
  coefficients = matrix(coefficients_vector, ncol = 3, byrow = TRUE)  #column = 3 since it is a "quadratic" function (degree of 2)
  
  #FUNCTIONS PER INTERVAL
  qsi.fxns = list()
  for (i in 1:no_of_intervals) {
    second_degree = paste(round(coefficients[i, 1], digits = 4), "x^2", sep = " * ")
    first_degree = paste(round(coefficients[i, 2], digits = 4), "x", sep = " * ")
    polynomial = paste(paste(second_degree, first_degree, sep = " + "), round(coefficients[i, 3], digits = 4), sep = " + ")
    qsi.fxns[[i]] = eval(parse(text = paste("function (x)", polynomial)))
  }
  
  #approximate x value using the functions (search for appropriate function using the intervals)
  for (i in 1:length(x_values)) {
    if(x >= x_values[i] && x <= x_values[i + 1]) {
      y = qsi.fxns[[i]](x)
      break
    }
  }
  
  values = list(qsi.fxns = qsi.fxns, y = y)
  return (values)
}

#USER INPUT
independentVector = c(3.0, 4.5, 7.0, 9.0)
dependentVector = c(2.5, 1.0, 2.5, 0.5)

data = list(independentVector, dependentVector)

result = poly.qsi(data, 5)
result