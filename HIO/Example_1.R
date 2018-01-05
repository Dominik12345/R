# This script implements an not HIO 2d-system with one observable

####################
##### R-Setup ######
####################

library(deSolve)



####################
##### Problem ######
####################

# DEFINE SYSTEM --->

# Set System Parameters
system.parameters <- c('alpha1' = 1, 'alpha2' = 1/2, 'beta21' = 1, 'beta12' = 0, 
                       'gamma1' = 1 , 'gamma2' = 0)

# Set Initial Values and Time 
state.parameters <- c('x1' = 1.3, 'x2' = 1.5)
observation.times <- seq(from = 0 ,to = 10,length.out = 100)

# Define System Equations 
system.dynamics = matrix( c(-system.parameters['alpha1'], system.parameters['beta21'], 
                            system.parameters['beta12'], -system.parameters['alpha2']),
                            nrow = 2 ,byrow = TRUE)

# Define Observation Matrix
system.observation = matrix( c(system.parameters['gamma1'], system.parameters['gamma2']),
                            nrow = 1 ,byrow = TRUE)

# <--- DEFINE SYSTEM

# ---> Define Hidden Inputs
temp.w1 <- function(t) {
  temp.out <-0#system.parameters['beta21']/system.parameters['alpha2'] * ( exp(-system.parameters['alpha2']*t) - 1)
  return(temp.out)
}
temp.w2 <- function(t) {
  temp.out <-0#1
  return(temp.out)
}


hi.input = c(temp.w1,temp.w2)
# <--- Define Hidden Inputs

# ---> Solve numerically
# Convert Into Standard Form
ODE.function <- function(t, state, parameters) {
  with(as.list(c(state,parameters)),{
    return(list(  c( system.dynamics %*% c(x1,x2))+c(hi.input[[1]](t),hi.input[[2]](t))  ))
    })
} 

ODE.solution <- ode( y = state.parameters, times = observation.times, func = ODE.function, 
                     parms = system.parameters)

observation.data <- data.frame('times' = observation.times)
observation.data['y'] <- 0
for ( i in seq( from = 1, to = length(ODE.solution[,1]), by =1) ) 
  observation.data$y[i] <- (system.observation %*% ODE.solution[i,-1])[[1]]

# <--- Solve numerically


####################
### Show Results ###
####################

plot(observation.data)