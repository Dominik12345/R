library(deSolve)

#define constants
kAlpha1 <- 1.5
kAlpha2 <- 0.5
kAlpha3 <- 0.3
kBeta1 <- 0.8
kBeta2 <- 0.2
kAmplitude <- 3
kOmega <- 1

x1.initial <- 5
x2.initial <-3
x3.initial <- 0

kTInitial <- 0
kTFinal <- 5
t.plot <- seq(kTInitial, kTFinal, by = 0.01)


# define system equations
f1 <- function(x1, x2, x3, u, w1 ,w2, w3){
  return( -kAlpha1 * x1 + u + w1)
}
f2 <- function(x1, x2, x3, u, w1, w2, w3){
  return( -kAlpha2 * x2 + kBeta1 * x1 + w2)
}
f3 <- function(x1, x2, x3, u, w1, w2, w3){
  return( -kAlpha3 * x3 + kBeta2 * x1 + w3)
}

#define known and hidden inputs
input.known <- function(t){
  return(t)
}
input.hidden.1 <- function(t){
  return( kAmplitude * sin(kOmega * t) )# t**2)
}
input.hidden.2 <- function(t){
  return( 0 )
}
input.hidden.3 <- function(t){
  return( 0 )
}

#define integrated functions
v1 <- function(t){
  return((input.known(t) + input.hidden.1(t)) * exp(kAlpha1 * t) )
}

V1 <- function(t){
  return(integrate( v1 , kTInitial, t)$val + x1.initial)
}

v2 <- function(t){
  return( ( input.hidden.2(t) * exp(kAlpha1 * t) + kBeta * sapply(t,V1)) * exp((kAlpha2-kAlpha1) * t) )
}

V2 <- function(t){
  return(integrate( v2 , kTInitial, t)$val + x2.initial)
}


#Exact Solutions 
x1 <- function(t){
  return(V1(t) * exp( - kAlpha1 * t))
}

x2 <- function(t){
  return(V2(t) * exp( - kAlpha2 * t))
}




#Solve Differential Equation Numerically
parameters <- c() #f1,f2,f3 are written as non-parametric model
state <- c(x1 = x1.initial, x2 = x2.initial, x3 = x3.initial)

f <- function(x1,x2,x3,t){
  output1.temp <- f1(x1,x2,x3,input.known(t), input.hidden.1(t), input.hidden.2(t),input.hidden.3(t))
  output2.temp <- f2(x1,x2,x3,input.known(t), input.hidden.1(t), input.hidden.2(t),input.hidden.3(t))
  output3.temp <- f3(x1,x2,x3,input.known(t), input.hidden.1(t), input.hidden.2(t),input.hidden.3(t))
  return(c(output1.temp,output2.temp, output3.temp))
}

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dx1 <- f(x1,x2,x3,t)[1]
    dx2 <- f(x1,x2,x3,t)[2]
    dx3 <- f(x1,x2,x3,t)[3]
    list(c(dx1,dx2,dx3))
  })
}

out <- ode(y = state, times = t.plot, func = Lorenz, parms = parameters)
head(out)

#PLOT Exact Solution
#dev.new()
par(mfrow = c(3,3)) #two plots in one window

'
x1.plot <- sapply(out[,1], x1)
x2.plot <- sapply(out[,1], x2)

plot(out[,1], x1.plot, xlab = "t", ylab ="",type = "n")
lines(out[,1], x1.plot, lty = 1)
lines(out[,1], x2.plot, lty =2)
legend("topleft", c("x1", "x2"), lty = c(1,2))
title("Analytical Solution")
'
# Plot Numerical solution
head(out[,1])
head(out[,2])
head(out[,3])
plot(out[,1], out[,2], xlab = "t", ylab = "" , lty = 1)
title("x1")
plot(out[,1], out[,3],  xlab = "t", ylab = "" ,lty = 1)
title("x2")
plot(out[,1], out[,4],  xlab = "t", ylab = "" ,lty = 1)
title("x3")
#legend("topleft", c("x1", "x2", "x3"), lty = c(1,2,3))
#title("Numerical Solution")

#Plot relative errors
'
plot(out[,1], (x1.plot-out[,2])/ x1.plot, xlab = "t", ylab ="(x1_exact-x1_numerical)/x1_exact" , type = "n")
lines(out[,1], (x1.plot-out[,2])/x1.plot , lty = 1)
title("Error x1")

plot(out[,1], (x2.plot-out[,3])/x2.plot ,xlab = "t", ylab = "(x2_exact-x2_numerical)/x2_exact", type = "n")
lines(out[,1], (x2.plot-out[,3])/x2.plot , lty =1)
title("Error x2")

'
#plot hidden inputs
w1.plot <- sapply(out[,1], input.hidden.1 )
w2.plot <- sapply(out[,1], input.hidden.2 )
w3.plot <- sapply(out[,1], input.hidden.3 )

plot(out[,1], w1.plot, xlab = "t", ylab = "",type = "n")
lines(out[,1],w1.plot, lty = 1)
title("w1")

plot(out[,1], w2.plot, xlab = "t" , ylab ="" ,type = "n")
lines(out[,1],w2.plot, lty = 1)
title("w2")

plot(out[,1], w3.plot, xlab = "t" , ylab ="" ,type = "n")
lines(out[,1],w3.plot, lty = 1)
title("w3")

# export data
output <- data.frame("t" = out[,1],"x1" = out[,2], "x2" = out[,3], "x3" = out[,4])

write.table(output, "C:/Users/kahl/Documents/R/Data/data.txt", sep = "\t")



print("Finished")