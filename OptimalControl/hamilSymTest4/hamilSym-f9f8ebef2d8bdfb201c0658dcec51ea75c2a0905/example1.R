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
  
  y1 = offset + x[,1] 
  y2 = offset + x[,2]
  y3 = offset + x[,3] 
  
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


#' create measurements
times <- seq(from= 0, to= 10, length.out = 30)
source('stateHiddenInput.R')
error <- function(t) 3*sin(t)
errorInd <- c(1,1,0,0)
errorMatrix <- data.frame(times,matrix(rep(error(times),3),ncol=3))
errorMatrix <- data.frame(mapply(`*`,errorMatrix,errorInd))

errorInput <- list()
errorInput$optW <- c(1,1,1)
errorInput$w <- apply(X = errorMatrix[,-1], MARGIN = 2, FUN = function(x) approxfun(x = errorMatrix[,1], y = x, method = 'linear', rule=2))
sol <- ode(y = x0, times = times, func = hiddenInputState, parms = parameters, input = errorInput)


X <- sol[,-1]
meas <- as.data.frame(do.call(cbind,testMessure(X)))
meas = as.data.frame(cbind(times,meas))

y <- meas

stepAlpha <- 4