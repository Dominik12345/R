# clear the Global Environment
rm(list=ls())


#############################
# 2.1 Vectors               #
#############################

# vectors store several values of the same mode

x <- c(1,2,3,4)
y <- 1:4
z <- rep(1,10)

words <- c('This', 'is', 'a', 'sentence', '.')

# What happens if you do ... ?

mixed.vector1 <- c(TRUE, 3.1)
mixed.vector2 <- c(TRUE, 'word', 3.7)

# What happens if you concatenate vectors ... ?

y <- c(mixed.vector1, mixed.vector2)

# Try to do basic operations and comparisons with vectors.
x <- c(2,4)
y <- c(1,2)

# Question: What is the difference between & and &&, | and || ?

# Note the recycling.

x <- 4:15
y <- c(1,-1)

x.alternating <- x * y

# We can grap the k-th component of a vector using [k]

x.alternating[2]

# or several components.

x.alternating[c(2,4,6)]

# Try to set only the first component of x.alternating to 'word'

# We get the number of components of vector x using length(x).

length(x.alternating)

# Question: What is the difference between 
#           x <- c() and x <- numeric(0) ?


#############################
# 2.2 Exercise              #
#############################

data <- rnorm(n = 1000, mean = 5, sd = 2)
plot(data)

# Problem:

# 1) Compute mean and standard deviation (std) 
#    of data 
# 2) grap the values in data that lie 
#    within mean +/- std and compute 
#    their relative number number


#############################
# 2.3 named vectors         #
#############################

values <- c(mean, std, relative.number)
names <- c('mean', 'std', 'rel num')

results <- values
names(results) = names

# Check results[1] and results['mean']

#############################
# 2.4 matrices              #
#############################

# Matrices are the two dimensional extension of vectors

A <- matrix( rep(TRUE,20), nrow = 5, ncol = 4)

B <- matrix(c(1,'word', TRUE, 3.1), nrow = 2)

C <- c(TRUE,FALSE)

# Try to grap the 1-2-component of B

# or the 1-st row/column of B.

#Note 
is.vector(B[1,])

# Question: What is the result of C * A ? ( A * C ?)

# We can combine matrices using cbind and rbind.
# Try it out.

A <- matrix( c('a','b','c','d'),nrow = 2)

# We can also perfrom matrix multiplication as known
# from linear algebra and invert square matrices

A <- matrix( c(0,1/2,1,0), nrow = 2)

# Question: What is the difference between A * A and A %*% A ?

# Try solve(A)

#############################
# 2.5 Exercise              #
#############################

# Consider the situation: 
# Three internet providers 'A' 'B' and 
# 'C' compete for customers.
# In January 'A' has got 10000 customers,
# 'B' 4000 and C '500'.
# Assume that after each month, 10% of the customers 
# of 'B' and 'C' change to 'A' and 5% of 'A' change 
# equally distributed to 'B' and 'C'.

# Problem:

# Compute the customer distribution in the next 3 years.
# Plot the number of customers of as a time series


#############################
# 2.6 lists                 #
#############################

# Lists can store several types of variables.

A <- matrix( rnorm(n = 6), nrow = 2)
names <- c('x1', 'x2', 'x3')
value <- 3

data <- list( 'matrix' = A, 'names' = names, 'value' = value, 'object' = NULL)

# We get items of list list using [[]] and $.
# Check it out.

#############################
# 2.7 Exercise              #
#############################

data <- list( 1:5 , TRUE , matrix( rnorm(n = 4), nrow = 2) , NULL, 
              list( 'values' = c(1,2,3)) , c(1,2,3,4,5) , 'test')


# Problem: 

# Label the elements of data with their mode, i.e.
# integer, numeric, logical, NULL, complex, character
# Write a short algorithm, that works for arbitrary lists.


#############################
# 2.8 data frames           #
#############################

# data.frames are lists of vectors.

data <- mtcars
print(data)

# they combine list structure
data$mpg
data[[1]]

# and matrix structure
data[6,1]
data['Valiant','mpg']

# Question: What is the difference between data[1] and data[[1]] ?

#############################
# 2.9 Exercise              #
#############################

data <- mtcars
val <-  rep('A',32)
'B' -> val[18] -> val[26]
data <- cbind.data.frame(data, val)

# Use class() and typeof() to get the classes and types of data and 
# of the columns of data.

# Use class conversion, e.g.

x <- 2.1
print( class(x) )
class(x) <- 'integer'
x <- as.factor(x)
x <- as.numeric(x)

# to convert the categorical variables of data into factors.

# Check out the classes and types again.

# Question: Explain the type of data$val. 