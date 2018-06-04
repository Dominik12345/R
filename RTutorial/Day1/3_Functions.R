###############################################################################
#                                                                             #
# R is made for 'functional programming'.                                     #
# That means                                                                  #
# - the objects interest, e.g. data is understood as 'object' of a            #
#   certain class,                                                            #
# - the actions and analysis are encapsulated within 'functions' and          #
#   'methods'.                                                                #
#                                                                             #  
###############################################################################

# clear the Global Environment
rm(list=ls())

#############################
# 3.1 Functions             #
#############################

# Functions and methods are predefined procedures.

# We define a function with a name 'TestFunction'
# that gets one 'argument' called 'word' within the 
# definition.

TestFunction <- function(word) {
  print(word)
  return(NULL)
}

# A 'function call' activates a function.
output <- TestFunction(word = 'Hello World!')


#############################
# 3.2 Exercise              #
#############################

data <- list()

for (i in 1:5) {
  data[[paste('y', as.character(i), sep = '')]] <- matrix(
          data = rnorm(n = 8**2, mean = i, sd = 0.1 * i), ncol = 8)
}

# Problem: Write a function, that accepts a list of same size numeric matrices and 
#          returns a dataframe showing the mean and standard deviation of 
#          the k-th column of each matrix.
#          Try to use only the already known basic operations and functions and 
#         nrow(A) to get the .

#############################
# 3.3 apply, sapply, lapply #
#############################

# Take the data 

particiant.name <- c('P1', 'P2', 'P3', 'P4')
gender <- c('m','f','f','m')
val1 <- rnorm(n = 4, mean = 1, sd = 0.2)
val2 <- rnorm(n = 4, mean = 10, sd = 5)

# and write it into a matrix and a data.frame.

# The functions apply() and sapply() simplify calculations that 
# are reapeated over matrices and data.frames.

# For matrices (and arrays) use apply and
# for data.frames use sapply and lapply.

# Questions: What is the difference between sapply and lapply ?

#############################
# 3.4 Exercise              #
#############################

data <- matrix( rnorm(n = 30000, mean = 1, sd = 0.2) ,
                ncol = 3)

# Take a look at data.

# Problem: Choose the maximum of x,y and z in each row.
#          Solve this 1) without apply 2) using apply and check 
#          if 1) and 2) really provide the same results.
#          Measure the computation time using system.time().
#          You can use the build-in functions max() and nrow() as 
#          well as all().

