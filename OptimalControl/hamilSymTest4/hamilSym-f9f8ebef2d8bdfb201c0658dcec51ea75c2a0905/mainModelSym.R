rm(list=ls())
graphics.off()

source("auxFunctions.R")
#setWd()

usePackage('Deriv')     # symbolic differentiation
usePackage('deSolve')   # solver for odes
usePackage('pracma')    # Practical Numerical Math Routines
usePackage("Hmisc")

source("classOdeEquation.R")
source("getEquations.R")
source("symbolicDiff.R")
source('errorBarSeg.R')
source('createMeassureData.R')

## generating test data
#source('example1.R')
#source('uvb.R')
#source('epo.R') 
source('Test_breit-wigner-resonance.R')


alpha1 <- 0.
alpha2 <- 0.001
#alpha2 <- 0.0001
stepAlpha <- 40

source('greedyApproach.R')
beta <- 0.8


results <- greedyApproach(alphaStep = stepAlpha,Beta = beta,alpha1 = alpha1, alpha2 = alpha2, 
                          x0 = x0, optW = optW , times=times, 
                          measFunc= testMessure,  measData = y, std = NULL,
                          parameters = parameters, 
                          modelFunc = testModel, greedyBoolean = TRUE)

#' Greedy Approach Algorithm
#' calculates controls based on a first optimisation with gradient descent; should result in a sparse vector
#' of hidden inputs
#' 
#' Arguments
#' alphaStep      the starting stepsize for the gradient descent
#'                a fitting stepsize will be calculated based on a backtracking line search
#'                if the algorithm converges to slow use a bigger stepsize
#' 
#' alpha1         L1-norm parameter of the dynamic elastic net approach, is set to zero for this algorithm
#' 
#' alpha2         L2-norm parameter of the dynamic elastic net approach
#'                used for regulation purposes
#'                set to NULL for a approximation of alpha2 - will results in a longer runtime
#'                
#' x0             inital state of the ODE system
#' 
#' optW           a vector that indicates for which knots of the network a input should be calculated
#' 
#' times          time sequence for which output is wanted; the first value of times must be the initial time
#' 
#' measFunc       a R-Function that is used for measurement of the states if the system is not completly 
#'                measurable; an empty argument will result in the assumption that the complete system is
#'                measurable
#'                
#' measData       a table that containts the measurements of the experiment; used to calculate the needed inputs
#' 
#' parameter      vector or named vector that contains the parameters of the ODE equation
#' 
#' modelFunc      a R-Function that states the ODE system for which the hidden inputs should be calculated
#' 
#' greedy         a boolean that states if the greedy approach should be used; if set to FALSE the algorithm 
#'                will only use perform a calculation of the inputs for all knots without a sparse solution
#'                
#' 
