parameters = 10.^c(-1.8,-0.71,-1.1,-0.54,-1.3,-3.2,-2.9,-1.9);
x0 = c(10^-0.54,10^-0.11,0,0,0,0)

testModel <- function(t,x,parameters) {
    k_t   = parameters[1]
    k_on  = parameters[2]
    k_off = parameters[3]
    b_max = parameters[4]
    k_e   = parameters[5]
    k_ex  = parameters[6]
    k_di  = parameters[7]
    k_de  = parameters[8]
  
      dx1 = k_t*b_max  - k_t*x[1] - k_on*x[1]*x[2]+k_off*x[3]+k_ex*x[4]
      dx2 = -k_on*x[1]*x[2]+k_off*x[3]+k_ex*x[4]
      dx3 = k_on*x[1]*x[2]-k_off*x[3]- k_e*x[3]
      dx4 = k_e*x[3]-k_ex*x[4]-k_di*x[4]-k_de*x[4]
      dx5 = k_di*x[4]
      dx6 = k_de*x[4]
    
    return(list(c(dx1,dx2,dx3,dx4,dx5,dx6)))
}


testMessure <- function(x) {
  
  offset = 1.0e-5
  scale = 1
  
  y1 = offset + x[,2] + x[,6]
  y2 = offset + x[,3]
  y3 = offset + x[,4] + x[,5]

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
times <- seq(0,300, length = 20)
times <- c(0, 10, 30, 60, 120, 180, 240, 300)
# generate input for testing that is reconstructed with dyn elastic net
error <- function(t) (1-1/(1+t))
errorInd <- c(1,0,0,0,0.1,0,0)
#errorInd <- c(1,0,0,0,0,0,0)
errorMatrix <- data.frame(times,matrix(rep(error(times),6),ncol=6))
errorMatrix <- data.frame(mapply(`*`,errorMatrix,errorInd))

errorInput <- list()
errorInput$optW <- c(1,1,1,1,1,1)
errorInput$w <- apply(X = errorMatrix[,-1], MARGIN = 2, FUN = function(x) approxfun(x = errorMatrix[,1], y = x, method = 'linear', rule=2))
source('stateHiddenInput.R')
sol <- ode(y = x0, times = times, func = hiddenInputState, parms = parameters, input = errorInput)


X <- sol[,-1]
meas <- as.data.frame(do.call(cbind,testMessure(X)))
meas = as.data.frame(cbind(times,meas))

y <- meas

stepAlpha = 0.0001