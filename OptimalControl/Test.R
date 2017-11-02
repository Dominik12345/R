#define constants

k1 <- 5 
k2 <- 5 
k3 <- 5 

a <- 1
b <- 1
c <- 1
d <- 10
e <- 0.1
eta <- 1

kTInitial <- 0
kTFinal <- 1
kNumberDatapoints <- 100
kDimOfMeasurement <- 1
kDimOfStatespace <- 2


# define functions ---

f.True <- function(x, u){
  output <- c( a*x[2] , c*x[3]+d*u , e*u )
  return( output )  
}

f.Nominal <- function(x, u){
  output <- c( a*x[2] , d*u )
  return( output )  
}

u <- function(t){
  return(exp(b*t))
}

x.True <- function(t){
  x1 <- a*c*e* 1/(b**3) * (u(t)-1) + a*d*1/(b**2) * (u(t)-1) +a*c*k3* 1/2*t**2 + k2*t+k1
  x2 <- c*e* 1/(b**2) * (u(t)-1) + d*1/(b) * (u(t)-1) +c*k3*t+ k2
  x3 <- e* 1/b * (u(t)-1) + k3
  return( c(x1,x2,x3) )
}

x.Nominal <- function(t){
  x1 <- a*d*1/b**2 * (u(t)-1)+t*k2+k1
  x2 <- d* 1/b * (u(t)-1) + k2
  return(c(x1,x2))
}

h <- function(x){
  output <- matrix(0, kDimOfMeasurement, 1)
  
  if (length(x) == 3) {
    output[,1] <- (eta * x[1])
  } else if (length(x) == 2) {
    output[,1] <- (eta * x[1])
  } else { 
    print("h: Problem with dimension of x")
    output[,1] <- (0)
  }
  return(output)
}


# --- define functions


# sample measurement ---

t.measured <- sort(runif(kNumberDatapoints, min = kTInitial, max = kTFinal), decreasing = FALSE)
x.measured <- sapply(t.measured, x.True)
y.measured <- matrix(0,1,kNumberDatapoints)

for (i in 1:kNumberDatapoints) {
  y.measured[1,i] <- h(x.measured[,i])
}

# --- sample measurement

# compute nominal-model expectations ---

x.nominal <- sapply(t.measured, x.Nominal)
y.nominal <- matrix(0,1,kNumberDatapoints)
for (i in 1:kNumberDatapoints) {
  y.nominal[1,i] <- h(x.nominal[,i])
}

t.plot <- seq(kTInitial, kTFinal, by = (kTFinal-kTInitial)/1000)
x.nominal.plot <- sapply(t.plot, x.Nominal)
y.nominal.plot <- matrix(0,1,length(t.plot))
for (i in 1:length(t.plot)) {
  y.nominal.plot[,i] <- h(x.nominal.plot[,i])
}

# --- compute nominal-model expectations


# plot measurent and nominal-expectation ---

plot(t.measured, y.measured[1,], xlab = "time / s", ylab = "y" )
lines(t.plot, y.nominal.plot[1,], lty = 1 )
title(main = "Measurement vs. Nominal Model")
legend('topleft',c("measured values", "nominal-model values"),pch = c(1,NA), lty = c(NA,1))

# --- plot measurement and nominal-expectation


# optimization ---

#weighting matrix and 
Q <- diag(kDimOfMeasurement)
kAlpha1 <- 1
kAlpha2 <- 1

# expected model and cost functional

f.Expected <- function(x, u, w) {
  return(f.Nominal(x, u) + w ) 
}

J <- function(w, y.measured, x, Q = diag(kDimOfMeasurement)){
  J <- 0
  for (k in 1:kNumberDatapoints) {
    for (i in 1:kDimOfMeasurement) {
      for (j in 1:kDimOfMeasurement) {
        J = J + Q[i,j] * (y.measured[i,k] - h(x[,k])[i,] ) * (y.measured[j,k]- h(x[,k])[i,] )
      }
    }
  
    for (i in 1:kDimOfStatespace) {
      J = J + kAlpha1 * sqrt(w[i,k]**2) + kAlpha2/2 * w[i,k]**2
    }
  }
  return(J)
}


# --- optimization

# TESTING AREA ---

# good guess of the hidden inputs
w.test <- matrix(0, kDimOfStatespace, kNumberDatapoints)
w.test[1,] <- x.measured[1,] - sapply(t.measured, x.Nominal)[1,] 
w.test[2,] <- x.measured[2,] - sapply(t.measured, x.Nominal)[2,] 

cost.test <- J(w.test,y.measured, x.nominal)
print(cost.test)
# --- TESTING AREA