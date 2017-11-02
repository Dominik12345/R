# usb network example

parameters = c(  ks1=0.23,
                 ks2=4.0526,
                 kdr1=0.1,
                 kdr2=0.2118,
                 k1=0.0043,
                 k2=161.62,
                 ka1=0.0372,
                 ka2=0.0611,
                 ka3=4.7207,
                 kd1=94.3524,
                 kd2=50.6973,
                 kd3=0.5508,
                 ks3=0.4397,
                 kdr3=1.246,
                 uv=1,
                 ka4=10.1285,
                 kd4=1.1999,
                 n1=3,
                 n2=2,
                 n3=3.5,
                 kdr3a=0.9735,
                 kdr3b=0.406,
                 ksr=0.7537,
                 FHY3_s=5)

x0 = c(0.2,10,2,0,0,20,0,0,0,4.2,0.25,20,0)

testModel <- function(t,x,parameters) {
  ks1 = parameters[1]
  ks2 = parameters[2]
  kdr1 = parameters[3]
  kdr2 = parameters[4]
  k1 = parameters[5]
  k2 = parameters[6]
  ka1 = parameters[7]
  ka2 = parameters[8]
  ka3 = parameters[9]
  kd1 = parameters[10]
  kd2 = parameters[11]
  kd3 = parameters[12]
  ks3 = parameters[13]
  kdr3 = parameters[14]
  uv = parameters[15]
  ka4 = parameters[16]
  kd4 = parameters[17]
  n1 = parameters[18]
  n2 = parameters[19]
  n3 = parameters[20]
  kdr3a = parameters[21]
  kdr3b = parameters[22]
  ksr = parameters[23]
  fhy3_s = parameters[24]
  
  dx1 = ((-2) * ((ka1 * (x[1]^2) * (x[4]^2)) - (kd1 * x[5])) + (-2) * ((ka2 * (x[1]^2) * x[2]) - (kd2 * x[3])) + ((ks1 *((1) + (uv * n3 * (x[11] + fhy3_s))))  - (kdr1 * ((1) + (n1 * uv)) * x[1])))
  dx2 = ((-1) * ((ka2*(x[1]^2) * x[2]) - (kd2 * x[3])) +(-1) * ((ka4 * x[2] * x[12]) - (kd4 * x[13])))
  dx3 = (((ka2 * (x[1]^2) * x[2]) - (kd2*  x[3]))) 
  dx4 = ((-2) * (k1*(x[4]^2)) + (2) * (k2 * x[6]) + (-2) * ((ka1 * (x[1]^2)* (x[4]^2)) - (kd1 * x[5])) + (-1)* (ka3 * x[4] *x[7]) + (1) * (2) * (k2 * x[6]))
  dx5 =  (((ka1 * (x[1]^2) * (x[4]^2)) -(kd1 * x[5])))
  dx6 = ((-1) * (k2 * x[6]) +  (k1 * (x[4]^2)) +(kd3 * (x[8]^2)))
  dx7 = ((-1) * (ka3 * x[4] * x[7]) + ((ks2 * ((1) + (uv * x[5]))) -(kdr2 * x[7])) + (2) * (kd3 * (x[8]^2)))
  dx8 = ((-2) * (kd3 * x[8]^2) + (ka3 * x[4] * x[7])) 
  dx9  = 0 
  dx10 = 0
  dx11 =  (((ks3 * ((1) + (n2 * uv))) -(kdr3 * (((x[3] / (kdr3a + x[3])) + (x[13] / (kdr3b + x[13]))) -(x[5] / (ksr + x[5]))) *  x[11])))
  dx12 = ((-1) * (ka4 * x[2] * x[12]) + (kd4 * x[13]))
  dx13 =((ka4 * x[2] * x[12]) - (kd4 * x[13]))
  
  return(list(c(dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx10,dx11,dx12,dx13)))
}


testMessure <- function(x) {
  
  y1 = 2*x[,5] + x[,4] + x[,8]
  y2 = 2*x[,5] + 2* x[,3] + x[,1]
  y3 = x[,6]
  y4 = x[,11]
  y5 = x[,4]
  
  return(list(y1,y2,y3,y4,y5))
}

source("classOdeEquation.R")
source("getEquations.R")
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
times <- seq(0, 0.25, length = 20)
times <- c(0, 0.01, .025, .05, 0.1,  0.15,  0.2, 0.25)
# generate input for testing that is reconstructed with dyn elastic net
error <- function(t) (1-1/(1+t))
errorInd <- c(1,0,0,1000,0,0,0,0,0,0,0,1000,0,0)
#errorInd <- c(1,0,0,0,0,0,0)
errorMatrix <- data.frame(times,matrix(rep(error(times),length(errorInd)-1),ncol=length(errorInd)-1))
errorMatrix <- data.frame(mapply(`*`,errorMatrix,errorInd))

errorInput <- list()
errorInput$optW <- rep(1,length(errorInd)-1)
errorInput$w <- apply(X = errorMatrix[,-1], MARGIN = 2, FUN = function(x) approxfun(x = errorMatrix[,1], y = x, method = 'linear', rule=2))
source('stateHiddenInput.R')
sol <- ode(y = x0, times = times, func = hiddenInputState, parms = parameters, input = errorInput)


X <- sol[,-1]
meas <- as.data.frame(do.call(cbind,testMessure(X)))
meas = as.data.frame(cbind(times,meas))

y <- meas

stepAlpha = 4
