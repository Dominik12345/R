rm(list=ls())
graphics.off()

source("auxFunctions.R")
usePackage('Deriv')
usePackage('deSolve')
usePackage('pracma')



source("classOdeEquation.R")
source("getEquations.R")
source("symbolicDiff.R")
source('errorBarSeg.R')
source('createMeassureData.R')

## generating test data
#source('example1.R')
#source('uvb.R')
#source('epo.R')
#source('C:/Users/kahl/Documents/R/OptimalControl/Controllability/Test_delta-peak.R')
#source('C:/Users/kahl/Documents/R/OptimalControl/Controllability/Test_linear.R')
#source('C:/Users/kahl/Documents/R/OptimalControl/Controllability/Test_linear2.R')
#source('C:/Users/kahl/Documents/R/OptimalControl/Controllability/Test_breit-wigner-resonance.R')
source('C:/Users/kahl/Documents/R/OptimalControl/Controllability/Test_breit-wigner-resonance_without_selfloop.R')

source('greedyApproach.R')
results <- greedyApproach(alphaStep = stepAlpha,x0 = x0, optW = rep(1,ncol(X)), times=times, measFunc= testMessure, measData = y,  parameters = parameters, modelFunc = testModel, odeObj = odeEq, greedyBoolean = TRUE)
