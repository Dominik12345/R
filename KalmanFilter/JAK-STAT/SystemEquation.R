# DEFINE GLOBAL PARAMETERS --->
# time (should be compared to data)
time.initial <- 0
time.final   <- 150
time.delta   <- 5
time.sequence<- seq(from = time.initial, 
                    to = time.final,
                    by = time.delta)

# system parameters from cite{Published 20 May 2010 on Science Express
# DOI: 10.1126/science.1184913}
system.parameters <- c(
  # in min^-1 and pM^-1
  "k_t"   = 0.0329366 , 
  "B_max" = 516 , 
  "k_on"  = 0.00010496 , 
  "k_off" = 0.0172135 , 
  "k_ex"  = 0.00993805 , 
  "k_e"   = 0.0748267 ,
  "k_di"  = 0.00317871 ,
  "k_de"  = 0.0164042 )
system.parameters <- append(system.parameters,c(
  # all equal 1 according to cite
  "psi_1" = 1 ,
  "psi_2" = 1 ,
  "psi_3" = 1 ))
system.parameters <- append(system.parameters,c(
  # initial values
  "x1.0" = system.parameters[["B_max"]],
  "x2.0" = 1,
  "x3.0" = 0.0,
  "x4.0" = 0.0,
  "x5.0" = 0.0,
  "x6.0" = 0.0 ))
system.parameters <- append(system.parameters,c(
  # relative std of the measurement noises
  "sigma1" = 0.1,
  "sigma2" = 0.1,
  "sigma3" = 0.1))


# <--- DEFINE GLOBAL PARAMETERS

# DEFINE SYSTEM --->
# dynamics "core model" from cite{Published 20 May 2010 on Science Express
# DOI: 10.1126/science.1184913}
system.Dynamics <- function(x,t,params,delta.t,...) {
  with(as.list(c(x,params)),{
    temp1 <- k_t*B_max - k_t * x1 - k_on*x1*x2 + k_off*x3 + k_ex*x4
    temp2 <- -k_on*x1*x2 + k_off*x3 + k_ex*x4
    temp3 <- k_on*x1*x2 - k_off*x3 - k_e*x3
    temp4 <- k_e*x3 - k_ex*x4 - k_di*x4 - k_de*x4
    temp5 <- k_di*x4
    temp6 <- k_de*x4
    out <- c("dx1" = temp1,
             "dx2" = temp2,
             "dx3" = temp3,
             "dx4" = temp4,
             "dx5" = temp5,
             "dx6" = temp6)
    return(out)
  })
}

# observation
system.ExpectedObservation <- function(x,t) {
  with(as.list(c(x , system.parameters)),{
    temp1 <- psi_1 * (x2 + 2 * x6)
    temp2 <- psi_2 * x3 
    temp3 <- psi_3 * (x4 + x5)
    out <- c( "y1" = temp1,
              "y2" = temp2,
              "y3" = temp3)
    return(out)
    })
}

system.ObservationCovariance <- function(params) {
  with(as.list(params),{
    temp <- matrix(c(sigma1,0,0,0,sigma2,0,0,0,sigma3),
                   nrow = 3,ncol = 3)
    return(temp)
  })
}

system.PerturbedObservation <- function(x,t,params,...) {
  with(as.list(c(x,params)),{
    temp0 <- system.ObservationCovariance(params = params)
    
    temp1 <- system.ExpectedObservation(x,t)
    
    temp11 <- temp1[["y1"]]
    temp21 <- temp1[["y2"]]
    temp31 <- temp1[["y3"]]

    temp12 <- rnorm(n = 1, mean = temp11, sd = temp0[1,1] * temp11)
    temp22 <- rnorm(n = 1, mean = temp21, sd = temp0[2,2] * temp21)
    temp32 <- rnorm(n = 1, mean = temp31, sd = temp0[3,3] * temp31)
    
    out <- c( "y1" = temp12,
              "y2" = temp22,
              "y3" = temp32)
    return(out)
  })
}

# <--- DEFINE SYSTEM


