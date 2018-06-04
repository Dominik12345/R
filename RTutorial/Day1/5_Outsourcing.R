###############################################################################
#                                                                             #
# R supports three object structures: S3, S4 and R5                           # 
#                                                                             #  
###############################################################################

# clear the Global Environment
rm(list=ls())

#############################
# 5.1 source files          #
#############################

source('filename')

# Exercise: Include the file 'myscript' that can be found in the folder 
#           'data'.

source('data/myscript.R')



#############################
# 5.2 libraries             #
#############################

# Install a library
install.packages('jsonlite')

# Load a library
library(jsonlite)

# Use funtions within the namespace '::' of a library
jsonlite::write_json( x = randu ,path = 'data/mydata')
