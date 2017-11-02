parameters = c(3, 4/(2*pi))
x0 = c(0,0,0,0)

testModel <- function(t,x,parameters) {
  a        = parameters[1];
  omega    = parameters[2];
  
  dx1 = a * sin(omega * t) * cos(omega * t)
  dx2 = x[1]
  dx3 = x[2]
  dx4 = x[3]
  
  return(list(c(dx1,dx2,dx3, dx4)))
}

testMessure <- function(x) {
  
  offset = 1.0e-5
  scale = 1
  
  y1 = x[,1] 
  y2 = x[,2]
  y3 = x[,3] 
  y4 = x[,4]
  
  return(list(y1,y2,y3,y4))
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


#' create measurements
times <- seq(from= 0, to= 10, length.out = 10)
source('stateHiddenInput.R')
error <- function(t) 10 * 1.1**2 / (1.1** 2 + (t - 5)**2)
errorInd <- c(0,0,1,0)
errorMatrix <- data.frame(times,matrix(rep(error(times),4),ncol=4))
errorMatrix <- data.frame(mapply(`*`,errorMatrix,errorInd))

errorInput <- list()
errorInput$optW <- c(1,1,1,1)
errorInput$w <- apply(X = errorMatrix[,-1], MARGIN = 2, FUN = function(x) approxfun(x = errorMatrix[,1], y = x, method = 'linear', rule=2))
sol <- ode(y = x0, times = times, func = hiddenInputState, parms = parameters, input = errorInput)


X <- sol[,-1]
meas <- as.data.frame(do.call(cbind,testMessure(X)))
meas = as.data.frame(cbind(times,meas))

y <- meas

stepAlpha <- 2