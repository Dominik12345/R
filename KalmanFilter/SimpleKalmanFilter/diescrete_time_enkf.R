# quick and dirty script to test a discrete time ensemble kalman filter

# LOAD PACKAGE ---> 
library(pomp)
# <--- LOAD PACKAGE 

# DEFINE THE DYNAMIC SYSTEM --->
time.initial <- 0
time.final <- 10
time.delta <- 0.1
time.sequence <- seq(from = time.initial, to = time.final,by = time.delta)
time.number <- length(time.sequence)

system.parameters <- c("x1.0" = 5 , "x2.0" = 5, "a" = 0.1, "b" = 0.1, 
                       "sigma1" = 0.5, "sigma2" = 0.5, "sigma3" = 0.5, "sigma4" = 0.5)

system.Dynamic <- function( x, t, params, delta.t,...) {
  with(as.list(c(x,params)), {
    temp1 <- -a * x1 * delta.t + x1   + rnorm(n = 1, mean = 0, sd = sigma1)
    temp2 <- (b * x1^2 - a * x2) * delta.t + x2   + rnorm(n = 1, mean = 0, sd = sigma2) 
    out <- c("x1" <- temp1, "x2" <- temp2)
    return(out)
  })
}
system.Observation <- function(x, t, params, ...) {
  with(as.list(c(x,params)),{
    temp1 <- x1 * x2 + rnorm(n = 1, mean = 0, sd = sigma3) 
    temp2 <- x1 + x2 + rnorm(n = 1, mean = 0, sd = sigma4)
    out <- c("y1" = temp1, "y2" = temp2)
    return(out)
  })
}

system.ExpectedObservation <- function(x,t) {
  with(as.list(x),{
    temp1 <- x1 * x2 
    temp2 <- x1 + x2 
    out <- c("y1" = temp1, "y2" = temp2)
    return(out)
  })
}

# <--- DEFINE THE DYNAMIC SYSTEM

# SAMPLE DATA --->

state.data <- data.frame("times" = time.sequence)
state.data["x1"] <- NA
state.data["x2"] <- NA
state.data[1,-1] <- c(system.parameters["x1.0"], system.parameters["x2.0"])
for (i in 2:length(time.sequence)) {
  state.data[i,-1] <- system.Dynamic(t = state.data$times[i-1],
                                     x = state.data[i-1,-1],
                                     params = system.parameters, 
                                     delta.t = time.delta )
}
observation.data <- data.frame("times" = state.data$times)
observation.data["y1"] <- NA
observation.data["y2"] <- NA
for (i in 1:length(time.sequence)) {
  observation.data[i,-1] <- system.Observation(x = state.data[i,-1],
                                               t = observation.data$times[i],
                                               params = system.parameters)
}
    


# <--- SAMPLE DATA
# INITIALIZE POMP OBJECT --->
pomp.initialize <- function(params, t0, ...) {
  with(as.list(params), {
    temp1 <- x1.0 + rnorm(n = 1, mean = 0, sd = sigma1)  
    temp2 <- x2.0 + rnorm(n = 1, mean = 0, sd = sigma2)
    out <- c("x1" = temp1, "x2" = temp2)
    return(out)
  })
}
pomp.object <- pomp( data = observation.data,
            times = "times", t0 = time.initial,
            rprocess = discrete.time.sim(step.fun = system.Dynamic, delta.t = time.delta),
            rmeasure = system.Observation,
            initializer = pomp.initialize,
            params = system.parameters
)

pomp.simulation <- simulate(pomp.object, nsim = 1)
# <--- INITIALIZE POMP OBJECT

# FILTER DATA --->
pomp.results <- enkf(object = pomp.object, #simulation, 
                     Np = 3,
                     verbose = FALSE,
                     params = coef(pomp.simulation),
                     h = system.ExpectedObservation,
                     R = matrix(c(coef(pomp.simulation)["sigma3"],0,0,
                                coef(pomp.simulation)["sigma4"]),ncol = 2)
                       )


# <--- FILTER DATA

# PLOT RESULTS --->
par(mfrow = c(2,2) )

plot(pomp.results@filter.mean[1,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[1,]    ,lty = 1)
lines(state.data$x1  ,lty = 2)
legend("topright", c("x1_filter", "x1_true"), lty = c(1,2))
title("x1")


plot(pomp.results@filter.mean[2,], xlab = "time step", ylab = "", 
     type = "n")
lines(pomp.results@filter.mean[2,]    ,lty = 1)
lines(state.data$x2  ,lty = 2)
legend("topright", c("x2_filter", "x2_true"), lty = c(1,2))
title("x2")


plot((pomp.results@filter.mean[1,]-state.data$x1)/state.data$x1, xlab = "time step", ylab = "",
     type = "n")
lines((pomp.results@filter.mean[1,]-state.data$x1)/state.data$x1, lty = 1 )
lines((pomp.results@filter.mean[2,]-state.data$x2)/state.data$x2, lty = 2 )
legend("bottomleft", c("Dx1", "Dx2"), lty = c(1,2))
title("relative errors")



plot( pomp.results@data[1,] , xlab = "time step", ylab = "",
     type = "n")
lines(pomp.results@data[1,] , lty = 1 )
lines(pomp.results@data[2,] , lty = 2 )
legend("topright", c("y1", "y2"), lty = c(1,2))
title("Observations")

# <--- PLOT RESULTS