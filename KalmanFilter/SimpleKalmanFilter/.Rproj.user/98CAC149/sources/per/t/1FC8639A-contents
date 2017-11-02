# some global parameters
time.initial <- 0 # initial time
time.final <- 20  # final time
time.delta <- 0.1
time.sequence <- seq(from = time.initial, to = time.final, by = time.delta) # measurement times

# DEFINE A CONTINUOUS TIME DYNAMIC SYSTEM --->

# set the initial state value 
system.parameters <- c("x1.0" = 1.10, "x2.0" = 2.5,
                       "a" = 0.1, "b" = 0.2,
                       "sigma1" = 0.1, "sigma2" = 0.1,
                       "sigma3" = 1.1, "sigma4" = 6.1)
# define the system dynamics with state x1,x2 and parameters a,b
system.Dynamic <- function( x, t, params, delta.t,...) {
  with(as.list(c(x,params)),{
  temp1 <-  x2   + rnorm(n = 1 , mean = 0, sd = sigma1)
  temp2 <- -x1 - a * x2    + rnorm(n = 1 , mean = 0, sd = sigma2)
  out <- c("x1" = temp1, "x2" = temp2)
  return(out)
  })
}

# define observation mapping and Covariance
system.ExpectedObservation<- function(x,t,params,...) {
  temp1 <- x[1]^2 + x[2]^2
  temp2 <- x[1]+x[2]
  out <- c("y1" = temp1, "y2" = temp2)
  return( out )
}

#error on observation about 10%
system.ObservationCovariance <- function(params) {
  with(as.list(params),{
    temp <- matrix(c(sigma3,0,0,sigma4),nrow = 2,ncol = 2)
    return(temp)
  })
  }

system.Observation<- function(x,t,params,...) {
  with(as.list(c(x,params)),{
  temp1 <- 100 + sqrt(x1^2 + x2^2) + rnorm(n = 1, mean = 0, sd = sigma3)
  temp2 <- 100 + x1+x2       + rnorm(n = 1, mean = 0, sd = sigma4)
  out <- c("y1" = temp1, "y2" = temp2)
  return( out )
  })
}

# <--- DEFINE A CONTINUOUS TIME DYNAMIC SYSTEM

pomp.initialize <- function(params, t0, ...) {
  with(as.list(params), {
    temp1 <- x1.0 + rnorm(n = 1, mean = 0, sd = sigma1)  
    temp2 <- x2.0 + rnorm(n = 1, mean = 0, sd = sigma2)
    out <- c("x1" = temp1, "x2" = temp2)
    return(out)
  })
}
