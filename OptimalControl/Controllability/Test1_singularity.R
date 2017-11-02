# This script generates data to test the greedy-approach dynamic elastic net
library(deSolve)


# The nominal system 
parameters = c(3, 4/(2*pi))
x0 = c(0,0,0,0)

testModel <- function(t,x,parameters) {
  a        = parameters[1];
  omega    = parameters[2];

  dx1 = a * sin(omega * t)**2
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




# Solve the true (perturbated) system with deSolve

# set the time data and the states to use it in deSolve
t.data <- seq(0, 10, length.out = 500)
parameters.solve <- c(a = parameters[1], omega = parameters[2], Gamma = 1.1, peak = 0, t.centre = mean(t.data))
state <- c( x1 = x0[1], x2 = x0[2] , x3 = x0[3], x4 = x0[4])

Perturbation <- function(t,x1,x2,x3,x4, Gamma, peak, t.centre) {
  #perturbation function
  perturbation.temp <- peak * Gamma**2 / (Gamma ** 2 + (t - t.centre)**2)
  
  return(perturbation.temp)
}

Lorenz <- function(t, state, parameters.solve) {
  with(as.list(c(state, parameters.solve)),{
    dx1 <- a * sin(omega * t) * cos(omega * t)
    dx2 <- x1
    dx3 <- x2 + Perturbation(t,x1,x2,x3,x4, Gamma, peak, t.centre)
    dx4 <- x3
    list(c(dx1,dx2,dx3, dx4))
  })
}

out <- ode(y = state, times = t.data, func = Lorenz, parms = parameters.solve)

# MAke plots

plot(out)

t.plot  <- out[,1]
x1.plot <- out[,2]
x2.plot <- out[,3]
x3.plot <- out[,4]
x4.plot <- out[,5]
y.plot <- rep(0, length(t.plot) )

for (i in 1:length(t.plot) ) {
  y.plot[i] <- Perturbation(t.plot[i], x1.plot[i], x2.plot[i], x3.plot[i], x4.plot[i], 
                            parameters.solve["Gamma"], parameters.solve["h"])
}

# how the perturbation should look like
#plot(t.plot, y.plot)

# export data
output <- data.frame("t" = out[,1],"x1" = out[,2], "x2" = out[,3], "x3" = out[,4], "x4" = out[,5])


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

source('stateHiddenInput.R')
stepAlpha <- 2

y <- output
X <- output
times <- output["t"]
times <- times[,1]
