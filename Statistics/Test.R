# Simulation of some Random Experiments

# Constants

kNMin <- 1
kNMax <- 10
knTotal <- 20
kOmega <- seq(kNMin, kNMax, by = 1)
kSetStrength <- 200

kp <- 0.7

# PROBABILITY FUNCTIONS --->
#Returns the probability of drawing x from Omega = {1,2,...,N}
GetProbabilityOfOneDiscreteEqual <- function(x, Omega = kOmega ){
  if ( !is.element(x,Omega) ) {
    print("Outcome is not an element of the sample space.")
    return(0)
  } else { return(1/length(Omega)) }
}


#Returns the probability of drawing x=(x1,x2,x3...,xn) from Omega = {1,2,...,N}
GetProbabilityOfDiscreteEqual <- function( x, Omega = kOmega){
  p <- 1
  for (i in x) {
    p.temp <- GetProbabilityOfOneDiscreteEqual(i,Omega)
    if (p.temp == 0) {
      print("Event is not a subset of the sample space.")
      return(0)
    } else {
      p <- p * p.temp
    }
  }
  return(p)
}

# <--- PROBABILITY FUNCTIONS

# RANDOM VARIABLES --->

MakeHistogram <- function(x, Omega = kOmega){
  hist.temp <- numeric( length = length(Omega) )
  for ( i in 1:length(x) ) {
    for (j in 1:length(Omega)) {
      if (x[i] == Omega[j]) {
        hist.temp[j] <- hist.temp[j] + 1
      }
    }
  }
  return(hist.temp)
}

BernoulliizeOne <- function(x, A){
  if (is.element(x,A)) {
    return(1)
  } else { return(0) }
}

Bernoulliize <- function(x, A) {
  bern.temp <- numeric(length(x))
    for (i in 1:length(x)) {
      bern.temp[i] <- BernoulliizeOne(x[i],A)
    }
  return(bern.temp)
}

Binomialize <- function(x, A) {
  bern.temp <- Bernoulliize(x,A)
  bin.temp <- MakeHistogram(bern.temp, c(0,1))
  return(bin.temp)
}

# <--- RANDOM VARIABLES


# GENERATE RANDOM NUMBERS --->

#equally distributet integers
x.equal <- sapply(runif(knTotal,kNMin,kNMax+1), floor)

#set of equally distributed integers
set.x.equal <- matrix(0, knTotal, kSetStrength)
for (i in 1:kSetStrength) {
  set.x.equal[,i] <- sapply(runif(knTotal,kNMin,kNMax+1), floor)
}

# binomial distributet numbers

A <- seq(kNMin, ceiling(kp * kNMax), by = 1)
binomial.k <- numeric(kSetStrength)
for (i in 1:kSetStrength) {
  binomial.k[i] <- Binomialize(set.x.equal[,i], A)[2]
}
binomial.hist <- MakeHistogram( binomial.k, seq( 0 , 1.2 * max(binomial.k), by = 1) ) 

# <--- GENERATE RANDOM NUMBERS



# TEST AREA



# PLOT AREA

#barplot(MakeHistogram(x.equal))
barplot(binomial.hist)
