# Kalman filtering algorithm
library(pomp)
observation.data <- read.table('observation.data.txt')
state.data <- read.table('state.data.txt')

# read model 
source("SystemEquation.R")

f <- function( x, t, params, delta.t,...) {
  with(as.list(c(x,params)),{
    temp <- system.Dynamic(x = x,t = t, params = params,delta.t = delta.t)
    temp = temp * delta.t + x
    out <- c("x1" = temp[[1]], "x2" = temp[[2]])
    return(out)
  })
}


# CONVERT DATA INTO POMP OBJECT --->
pomp.object <- pomp(
  data   = observation.data        ,
  times  = "time"                  ,
  t0     = observation.data$time[1],
  params = system.parameters,
  rprocess = euler.sim(step.fun = f, delta.t = time.delta),
  rmeasure = system.Observation,
  initializer = pomp.initialize,
)
# <--- CONVERT DATA INTO POMP OBJECT

# SAMPLE ANOTHER DATA SET USING POMP --->
#new.data <- simulate(po)

# <--- SAMPLE ANTOTHER DATA SET USING POMP

# APPLY ENSEMBLE KALMAN FILTER --->

pomp.results <- enkf(object = pomp.object,
                    Np = 36,
                    verbose = FALSE,
                    params = coef(pomp.object),
                    h = system.ExpectedObservation,
                    R = system.ObservationCovariance(params = coef(pomp.object))
                    )

# <--- APPLY ENSELMBLE KALMAN FILTER

# PLOT --->
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


ymax = max(pomp.results@data[1,],pomp.results@data[2,] )
ymin = min(pomp.results@data[1,],pomp.results@data[2,] )

plot( pomp.results@data[1,] , xlab = "time step", ylab = "",
      type = "n", ylim = c(ymin,ymax) )
lines(pomp.results@data[1,] , lty = 1 )
lines(pomp.results@data[2,] , lty = 2 )
legend("topright", c("y1", "y2"), lty = c(1,2))
title("Observations")

# <--- PLOT
