#Einfaches Beispiel zur lasso und ridge Regularisierung


# DEFINE SOME CONSTANTS

kDimensionOfStatespace <- 2
kDimensionOfObservationspace <-2
# time interval
kTInitial <- 0
kTFinal   <- 10
kNumberOfTimesteps <- 100
kNumberOfPreOptimizationsteps <- 20
# regularisation parameters
# 1 = lasso, 2 = ridge
kAlpha1 <- 2.2
kAlpha2 <- 2.0
# noise amplitude
kNoiseAmplitude <- 0.1
kNumberOfRandomsteps <- 1000



# DEFINE A TRUE SYSTEM -----------------------

# state evolution in time
State1True <- function(t) {
  return( 3 * t )
}
State2True <- function(t) {
  return( 5 * t**2 )
}

# define outputmatrix
observation.matrix <- matrix(0,  nrow = kDimensionOfStatespace, ncol = kDimensionOfObservationspace )

observation.matrix[1,1] <- 1
observation.matrix[2,2] <- 1

# define timesequence
time <- seq( kTInitial, kTFinal, length.out = kNumberOfTimesteps)

#----------------------------------------------

  
# DEFINE AN UNFIT SYSTEM ---------------------
# w is kDimensionOfStatespace x number of functions 
# define the set of possibe functions
VectorOfFunctions <- function(t) {
  return(c(t, t**2, t**3))
}

StateUnfit <- function(w,t) {
  return( w %*% VectorOfFunctions(t) )
}

ObservationUnfit <- function(w,t) {
  output <- matrix(0, ncol = kDimensionOfObservationspace, nrow = length(t))
  for (i in 1:length(t)) {
    output[i,] <- observation.matrix %*% StateUnfit(w,t[i]) 
  }  
  return(output)   
}
#----------------------------------------------

# CREATE MEASURED DATA -----------------------
noise.state1.true = runif( kNumberOfTimesteps, -kNoiseAmplitude, kNoiseAmplitude)
noise.state2.true = runif( kNumberOfTimesteps, -kNoiseAmplitude, kNoiseAmplitude)
# evolution of states in time
state1.true <- sapply(time, State1True) 
state2.true <- sapply(time, State2True) 

# compute observation = (time, y1, y2, ...)
observation <- matrix(0, ncol = kDimensionOfObservationspace +1 , nrow = kNumberOfTimesteps) 
for (i in 1:kNumberOfTimesteps) {
  observation[i,] <-c(time[i], observation.matrix %*% c(state1.true[i], state2.true[i]) )
}

colnames.temp <- c("t" ,"y1", "y2")
dataframe.true <- as.data.frame(observation)
colnames(dataframe.true) <- colnames.temp
#----------------------------------------------

# DEFINE COST FUNCTIONAL --------------------
Regularization <- function(w) {
  lasso <- 0
  ridge <- 0
  for (i in 1:length(w)) {
    lasso = lasso + abs(w[i])
    ridge = ridge + abs(w[i])**2
  }
  lasso = kAlpha1 * lasso
  ridge = 1./2. * kAlpha2 * ridge 
}

TotalCost <- function(w, data) {
  cost.temp <- 0
  for (i in 1:kNumberOfTimesteps) {
    cost.temp = cost.temp + (ObservationUnfit(w,data[,"t"])[i,1]-data[i,"y1"] )**2 + 
        (ObservationUnfit(w,data[,"t"])[i,2] - data[i,"y2"])**2 + Regularization(w)
  }
  return( cost.temp )
}
#----------------------------------------------

# Finde optimal w by randomly varying w
# this should work since Regularization(w) is convex
w.initial <- matrix(0, nrow = kDimensionOfStatespace, ncol = length(VectorOfFunctions(0) ) )

temp <- matrix(runif(6,0,1),3,2)

Optimize <- function(w.initial, data) {
  #initial cost
  cost.initial <- TotalCost(w.initial, data)
  
  #define Algorithm to find an appropriate w
  w.process <- w.initial
  cost.process <- cost.initial
  
    
  for (i in 1:kNumberOfRandomsteps) {
    dw = matrix(runif(6,-0.1,0.1), nrow = kDimensionOfStatespace, ncol = length(VectorOfFunctions(0) ) )
    if(TotalCost(w.process + dw, data) < cost.process) {
      w.process = w.process + dw
      cost.process = TotalCost(w.process,data)
    }
  }
  print("w Matrix")
  print(w.process)
  print("costs")
  print(cost.process)
  return(w.process)
}


w.initial <- matrix(c(3,1,3,5,0,2),2,3)

w.test <- Optimize(w.initial, dataframe.true)
