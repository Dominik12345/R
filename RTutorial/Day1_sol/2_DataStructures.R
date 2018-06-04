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

# R converts different data types automatically

mixed.vector1 <- c(TRUE, 3.1)
mixed.vector2 <- c(TRUE, 'word', 3.7)

# we can concatenate vectors

y <- c(mixed.vector1, mixed.vector2)

# and perform basic operations componentwise

x <- c(2,4)
y <- c(1,2)

x + y
x - y
x * y
x / y

x ** 2

x > 3

# note the recycling

x <- 4:15
y <- c(1,-1)

x.alternating <- x * y

# we can grap one k-th component of a vector using [k]

x.alternating[2]
x.alternating[c(2,4,6)]

x.alternating[1] <- 'word'

# we get the number of components of vector x using length(x)

length(x.alternating)

# Question: What is the difference between 
#           x <- c() and x <- numeric(0) ?

x <- c()
x <- numeric(0)

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

# Possible solution:

# 1)
# iterate over index
mean <- 0
for (i in (1:length(data)) ) {
  mean = mean + data[i] / length(data)
}

# or iterate over components
mean <- 0
for ( value in data ) {
  mean = mean + value / length(data)
}

std <- 0
for (value in data) {
  std = std + (value - mean)**2 / (length(data) - 1)
}
std = std**(1/2)

# 2)

data.logical <- data > mean - std & data < mean + std

data.subset <- data[data.logical]

plot(data.subset)

relative.number <- length(data.subset) / length(data)


#############################
# 2.3 named vectors         #
#############################

values <- c(mean, std, relative.number)
names <- c('mean', 'std', 'rel num')

results <- values
names(results) = names

results[1]
results['mean']



#############################
# 2.4 matrices              #
#############################

# matrices are the two dimensional extension of vectors

A <- matrix( rep(TRUE,20), nrow = 5, ncol = 4)

B <- matrix(c(1,'word', TRUE, 3.1), nrow = 2)

C <- c(TRUE,FALSE)

# we can get the i,j - component
B[1,2]

# or the i-th row/column
B[1,]
B[,1]

#note 
is.vector(B[1,])

# Question: What is the result of C * A ? ( A * C ?)

A * C

C * A

# i.e. in doubt, matrices are treated as vectors

# We can combine matrices
A <- matrix( c('a','b','c','d'),nrow = 2)

cbind(A,A)
rbind(A,A)

# We can also perfrom matrix multiplication as known
# from linear algebra and invert square matrices

A <- matrix( c(0,1/2,1,0), nrow = 2)

A * A   # componentwise matrix multiplication
A %*% A # algebraic matrix multiplication

solve(A) # inverse, if possible


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

# Possible solution:

customers <- c('A' = 10000, 'B' = 4000, 'C' = 500)
customers.A <- c(as.numeric( customers['A'] ))

M <- matrix( data = c(  0.95, 0.1, 0.1,
                       0.025, 0.9,   0,
                       0.025,   0, 0.9 ),
                       nrow = 3, byrow = TRUE)
rownames(M) <- c('A','B','C')
colnames(M) <- c('A','B','C')

for (i in 1:36) {
  customers = M %*% customers
  customers.A = c(timeseries, customers[1])
}
plot(customers.A)




#############################
# 2.6 lists                 #
#############################

# lists can store several types of variables

A <- matrix( rnorm(n = 6), nrow = 2)
names <- c('x1', 'x2', 'x3')
value <- 3

data <- list( 'matrix' = A, 'names' = names, 'value' = value, 'object' = NULL)

data[[1]]
data$matrix

#############################
# 2.7 Exercise              #
#############################

data <- list( 1:5 , TRUE , matrix( rnorm(n = 4), nrow = 2) , NULL, 
              list( 'values' = c(1,2,3)) , c(1,2,3,4,5) , 'test')


# Problem: 

# Label the elements of data with their mode, i.e.
# integer, numeric, logical, NULL, complex, character
# Write a short algorithm, that works for arbitrary lists.

# Possible solution:

for (i in 1:length(data)) {
  names(data)[i] = class(data[[i]])
}


#############################
# 2.8 data frames           #
#############################

# data frames are lists of vectors

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

print( class(data) )
for (label in names(data)) {
  print(paste(label,class(data[[label]]) , sep = ' is object of class '))
  
}

print( typeof(data) )
for (label in names(data)) {
  print(paste(label,typeof(data[[label]]) , sep = ' is object of type '))
  
}


# Use class conversion, e.g.

x <- 2.1
print( class(x) )
class(x) <- 'integer'
x <- as.factor(x)
x <- as.numeric(x)

# to convert the categorical variables of data into factors.

data$cyl <- as.factor(data$cyl)
data$vs <- as.factor(data$vs)
data$am <- as.factor(data$am)
data$gear <- as.factor(data$gear)
data$carb <- as.factor(data$carb)

# Check out the classes and types again.

print( class(data) )
print( typeof(data)) 
for (label in names(data)) {
  print(paste(label,class(data[[label]]) , sep = ' is object of class '))
  print(paste(label,typeof(data[[label]]) , sep = ' is object of type '))
}

# Question: Explain the type of data$val. 