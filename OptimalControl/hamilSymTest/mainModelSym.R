rm(list=ls())
graphics.off()

source("auxFunctions.R")
#setWd()
usePackage('Deriv')
usePackage('deSolve')
usePackage('pracma')



source("classOdeEquation.R")
source("getEquations.R")
source("symbolicDiff.R")
source('errorBarSeg.R')
source('createMeassureData.R')

## generating test data
parameters = c(1.5, 0.5, 0.3, 0.8, 0.2)
x0 = c( 5 , 3 ,0)

testModel <- function(t,x,parameters) {
  alpha1   = parameters[1];
  alpha2   = parameters[2];
  alpha3   = parameters[3];
  beta1    = parameters[4];
  beta2    = parameters[5];

  dx1 = -alpha1 * x[1] + t
  dx2 = -alpha2 * x[2] + beta1 * x[1]
  dx3 = -alpha3 * x[3] + beta2 * x[1]
  
  return(list(c(dx1,dx2,dx3)))
}

testMessure <- function(x) {
  
  offset = 1.0e-5
  scale = 1
  
  y1 =  x[,1] 
  y2 =  x[,2]
  y3 =  x[,3]
  
  return(list(y1,y2,y3))
}

## create the files and functions
odeEq <- odeEq()

# get the model equations
odeEq <- createModelEqClass(odeEq,testModel)

# calculate the costates and jacobian
odeEq <- setMeassureFunc(odeEq,testMessure)
odeEq <- isDynElaNet(odeEq)
source('symbolicDiff.R')
odeEq <- calculateCostate(odeEq)

odeEq
source("createFunctions.R")
createFunctions(odeEq)







#' generate noisy data for testing purposes
#' ##########################################################
## integrate and add noise to measurements

data <- read.table("C:/Users/kahl/Documents/R/Data/data.txt", header=TRUE,
                             sep="\t")
times <- data["t"]
#y <- subset(data, select = c("t","x1","x2"))
y <- subset(data, select = c("t", "x1", "x2" , "x3"))
names(y) <- NULL
names(times) <- NULL
times <- times[,1]


source("dynElasticNet.R")
source("stateHiddenInput.R")
results <- dynElasticNet(x0 = x0, times=times, measFunc= testMessure, measData = y,  parameters = parameters, modelFunc = testModel, odeObj = odeEq)
results$optW

#barplot(colSums(results$w))

par(mfrow = c(3,3)) 
plot(results$x[1:100],results$x[101:200]           ,xlab = "t", ylab = "") ; title("x1")
plot(results$x[1:100],results$x[201:300]           ,xlab = "t", ylab = "") ; title("x2")
plot(results$x[1:100],results$x[301:400]           ,xlab = "t", ylab = "") ; title("x3")
plot(results$x[1:100], results$w[1:100]  ,xlab = "t", ylab = "") ; title("w1") 
plot(results$x[1:100], results$w[101:200],xlab = "t", ylab = "") ; title("w2")
plot(results$x[1:100], results$w[201:300],xlab = "t", ylab = "") ; title("w3")
