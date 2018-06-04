###############################################################################
#                                                                             #
# Run current line: Ctrl + Enter / Strg + Enter                               #
# Run script: Ctrl + Alt + R / Strg + Alt + R                                 #
#                                                                             #
# "#" makes text invisible "comment"                                          #
#                                                                             #    
# For styleguide: https://google.github.io/styleguide/Rguide.xml              #                                        
#                                                                             #
#################### maximal line length = 80 chars ###########################

# clear the Global Environment
rm(list=ls())

#############################
# 1.1 Basic Data Types      #
#     and Exmaples          #      
#############################


# NULL

x <- NULL

# variable x is initialized with the value NULL

mode(x)

# logicals

A <- TRUE
B <- FALSE

A & B

A | B

mode(A)

# numerics

x <- 7

y <- x

y <- 3.9

mode(x)

# integers

is.integer(x)

x <- as.integer(x)
y <- as.integer(y)

# complex numbers

z <- 1 + 1i

mode(z)

# characters

word1 <- 'one'
word2 <- 'two'

word3 = paste(word1, word2)

sub('wo', 'hree', w)

# Question: What is the value of word3 ?


#############################
# 1.2 Basic Operations      #      
#############################

# numeric calculation

x <- 3
y <- 4

x + y ; x * y # ";" separates statements
x - y 
x / y


x*2 # skalar multiples

x**2 # powers

x^2 # powers

# comparisons return logicals

x <- 3
y <- 'Hans'

x == y

x != y

y <- 4

x > y
x >= y

#############################
# 1.3 If statements         #      
#############################

A <- TRUE
B <- FALSE

if (A) {
  print('Hello world!')
}

if (B) {
  print('yes')
} else {        # syntax !
  print('no')
}


#############################
# 1.4 For loop              #      
#############################


for (k in (1:10) ) {
  print(paste('step', k))  
}


#############################
# 1.5 Exercise              #  
#############################

# Problem:

# Implement the factorial of n
# i.e. n * (n-1) * ... * 2 * 1
# using if/for statements
# note: n has to be integer

n <- 7.25

# Possible solution:

if (!is.integer(n)) {
  print('n has to be integer')
} else {
  
  factorial <- 1
  for (k in (1:n)) {
    factorial = factorial * k
  }
  print(factorial)
}





