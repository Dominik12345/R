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

########
# NULL #
########

# Initialize a variable x with the value NULL.

# Check mode(x)

############
# logicals #
############

#Initialize A and B with the values TRUE and FALSE.

# Check mode(A)

# Which effect have the operatons A & B and A | B ?

############
# numerics #
############

# Initialize x with the value 7 and check mode(x)

# Initialize y with the value x

# Change the value if y to 3.1. 
# Did you notice something?

############
# integers #
############

# Check is.integer(x)

# Now convert x and y into integers using as.integer(x)

###################
# complex numbers #
###################

# Initialize z with the value 1 + 1i

# Check mode(z)

##############
# characters #
##############

# Initialize word1 and word2 with 'one' and 'two'

# Check mode(word1)

word3 = paste(word1, word2)

sub('wo', 'hree', w)

# Question: What is the value of word3 ?


#############################
# 1.2 Basic Operations      #      
#############################

# numeric calculation

x <- 3
y <- 4

# Try to do basic numeric calculations.

# Try to compute powers.

# Try to do logical comparisons.

x <- 3
y <- 'Hans'


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


