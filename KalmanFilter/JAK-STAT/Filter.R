# INCLUDE LIBRARY AND DATA --->
# load pomp library
library(pomp)
# load JSON library to export data
library(rjson)
# get data 
observation.data <- read.table('observation.data.txt')
unperturbedobservation.data <- read.table('unperturbedobservation.data.txt')
state.data <- read.table('state.data.txt')
#read model
source("SystemEquation.R")
# <--- INCLUDE LIBRARY AND DATA 


# DEFINE FILTER PARAMETERS --->
pomp.Covariance <- function(params) {
  with(as.list(params),{
    out <- matrix(c(
      sigma1 , 0 , 0 ,
      0 , sigma2 , 0 ,
      0 , 0 , sigma3),
      byrow = TRUE,
      ncol = 3, 
      nrow = 3)
  })
}

pomp.parameters <- c(
  # std of the random shift in the stochastic EnKF
  "sigma1" = 1,
  "sigma2" = 1,
  "sigma3" = 1)

pomp.Initialize <- function(params, t0, ...) {
  with(as.list(params),{
    temp1 <- x1.0 + rnorm( n = 1, mean = 0,  sd = 0.1 * x1.0) 
    temp2 <- x2.0 + rnorm( n = 1, mean = 0,  sd = 0.1 * x2.0)  
    temp3 <- x3.0 + rnorm( n = 1, mean = 0,  sd = 0.1 * x3.0)  
    temp4 <- x4.0 + rnorm( n = 1, mean = 0,  sd = 0.1 * x4.0)  
    temp5 <- x5.0 + rnorm( n = 1, mean = 0,  sd = 0.1 * x5.0)  
    temp6 <- x6.0 + rnorm( n = 1, mean = 0,  sd = 0.1 * x6.0) 
    out <- c( "x1" = temp1, 
              "x2" = temp2, 
              "x3" = temp3,  
              "x4" = temp4,   
              "x5" = temp5,    
              "x6" = temp6) 
    return(out) 
  })
}


pomp.Dynamics <- function(x,t,params,delta.t,...) {
  with(as.list(c(x,params)),{
    temp <- system.Dynamics(x = x,
                            t = t,
                            params = params,
                            delta.t = delta.t) 
    out <-  x + temp * delta.t
    return(out)
  })
}
# <--- DEFINE FILTER PARAMETERS


# CREATE POMP OBJECT AND APPLY ENKF --->
pomp.object <- pomp(
  data = observation.data ,
  times = "time" ,
  t0 = observation.data$time[1],
  params = system.parameters,
  #initializer = pomp.Initialize,
  rprocess = euler.sim(step.fun = pomp.Dynamics, delta.t = time.delta) ,
  rmeasure = system.PerturbedObservation
)

pomp.results <- enkf(
  object = pomp.object ,
  Np = 10 ,
  verbose = FALSE ,
  params = coef(object = pomp.object),
  h = system.ExpectedObservation,
  R = pomp.Covariance(params = pomp.parameters)
)



# <--- CREATE POMP OBJECT AND APPLY ENKF


# WRITE SAMPLED DATA INTO FILE --->
out <- t(as.matrix(pomp.results@filter.mean))
out <- as.data.frame(out)
out$time <- rownames(out)

write.table( out, 'state.filter.txt', sep = "\t")
write( toJSON(out)   , 
       file = "C:/Users/kahl/Documents/Data/review_filterstate_json.txt")
# <--- WRITE SAMPLED DATA INTO FILE


# PLOT RESULTS --->
pdf(file = "C:/Users/kahl/Documents/Images/review_R.pdf")

par(mfrow = c(3,3) )

plot(observation.data$y1, xlab = "time step", ylab = "", type = "n")
points(observation.data$y1, lty = 1)
lines(unperturbedobservation.data$y1, lty =2)
lines(pomp.results@forecast[1,], lty = 3)
legend("topright", c("y1", "y1 unpert", "y1 forecast"), lty = c(1,2,3))
title("y1")

plot(observation.data$y2, xlab = "time step", ylab = "", type = "n")
points(observation.data$y2, lty = 1)
lines(unperturbedobservation.data$y2, lty =2)
lines(pomp.results@forecast[2,], lty = 3)
legend("topright", c("y2", "y2 unpert", "y2 forecast"), lty = c(1,2,3))
title("y2")

plot(observation.data$y3, xlab = "time step", ylab = "", type = "n")
points(observation.data$y3, lty = 1)
lines(unperturbedobservation.data$y3, lty =2)
lines(pomp.results@forecast[3,], lty = 3)
legend("topright", c("y3", "y3 unpert", "y3 forecast"), lty = c(1,2,3))
title("y3")




plot(pomp.results@filter.mean[1,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[1,],lty = 1)
lines(state.data$x1  ,lty = 2)
legend("topright", c("x1_filter", "x1_true"), lty = c(1,2))
title("x1")

plot(pomp.results@filter.mean[2,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[2,],lty = 1)
lines(state.data$x2  ,lty = 2)
legend("topright", c("x2_filter", "x2_true"), lty = c(1,2))
title("x2")

plot(pomp.results@filter.mean[3,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[3,],lty = 1)
lines(state.data$x3  ,lty = 2)
legend("topright", c("x3_filter", "x3_true"), lty = c(1,2))
title("x3")

plot(pomp.results@filter.mean[4,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[4,],lty = 1)
lines(state.data$x4  ,lty = 2)
legend("topright", c("x4_filter", "x4_true"), lty = c(1,2))
title("x4")

plot(pomp.results@filter.mean[5,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[5,],lty = 1)
lines(state.data$x5  ,lty = 2)
legend("topright", c("x5_filter", "x5_true"), lty = c(1,2))
title("x5")

plot(pomp.results@filter.mean[6,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[6,],lty = 1)
lines(state.data$x6  ,lty = 2)
legend("topright", c("x6_filter", "x6_true"), lty = c(1,2))
title("x6")

dev.off()

# <--- PLOT RESULTS

# CALCULATE ACCURACY --->
errors.absolute <- data.frame(
  "Dx1" = pomp.results@filter.mean[1,] - state.data$x1 ,
  "Dx2" = pomp.results@filter.mean[2,] - state.data$x2 ,
  "Dx3" = pomp.results@filter.mean[3,] - state.data$x3 ,
  "Dx4" = pomp.results@filter.mean[4,] - state.data$x4 ,
  "Dx5" = pomp.results@filter.mean[5,] - state.data$x5 ,
  "Dx6" = pomp.results@filter.mean[6,] - state.data$x6 )
errors.relative <- errors.absolute[-1,] / state.data[-1,-1]
errors.absolute.mean <- c(
                 "x1" =  mean(errors.absolute$Dx1), 
                 "x2" =  mean(errors.absolute$Dx2), 
                 "x3" =  mean(errors.absolute$Dx3), 
                 "x4" =  mean(errors.absolute$Dx4), 
                 "x5" =  mean(errors.absolute$Dx5), 
                 "x6" =  mean(errors.absolute$Dx6)) 
errors.absolute.std  <- c(
                 "x1" =  sd(errors.absolute$Dx1), 
                 "x2" =  sd(errors.absolute$Dx2), 
                 "x3" =  sd(errors.absolute$Dx3), 
                 "x4" =  sd(errors.absolute$Dx4), 
                 "x5" =  sd(errors.absolute$Dx5), 
                 "x6" =  sd(errors.absolute$Dx6)) 

errors.relative.mean  <- c(
                 "x1" =  mean(errors.relative$Dx1), 
                 "x2" =  mean(errors.relative$Dx2), 
                 "x3" =  mean(errors.relative$Dx3), 
                 "x4" =  mean(errors.relative$Dx4), 
                 "x5" =  mean(errors.relative$Dx5), 
                 "x6" =  mean(errors.relative$Dx6)) 
errors.relative.std  <- c(
                 "x1" =  sd(errors.relative$Dx1), 
                 "x2" =  sd(errors.relative$Dx2), 
                 "x3" =  sd(errors.relative$Dx3), 
                 "x4" =  sd(errors.relative$Dx4), 
                 "x5" =  sd(errors.relative$Dx5), 
                 "x6" =  sd(errors.relative$Dx6)) 
errors.out <- matrix( c( errors.absolute.mean, errors.absolute.std,
                         errors.relative.mean, errors.relative.std),
                      ncol = 4 )
dimnames(errors.out) = list( c("x1", "x2", "x3", "x4", "x5", "x6") ,
                             c("Mean Error","Mean Standard Deviation",
                               "Mean Relative Error", "Mean Relative Standard Deviation")
                             )
print( t(errors.out))
# <--- CALCULATE ACCURACY