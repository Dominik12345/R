# LOAD LIBRARIES AND SOURCES --->
library(rjson)
library(deSolve)
source('SystemEquation.R')

# <--- LOAD LIBRARIES AND SOURCES

# SOLVE IVP USING deSolve --->

IVP.Dynamics <- function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    dx <- system.Dynamics(x = state,
                          t = t,
                          params = parameters,
                          delta.t = time.delta)
    return(list(dx))
    
  })
}
IVP.initialstate <- c( "x1" = system.parameters[["x1.0"]],
                       "x2" = system.parameters[["x2.0"]],
                       "x3" = system.parameters[["x3.0"]], 
                       "x4" = system.parameters[["x4.0"]], 
                       "x5" = system.parameters[["x5.0"]], 
                       "x6" = system.parameters[["x6.0"]])
IVP.solution <- ode( y = IVP.initialstate,
                     times = time.sequence,
                     func = IVP.Dynamics,
                     parms = system.parameters)
    
IVP.solution <-  data.frame("time" = IVP.solution[,1],
                            "x1" = IVP.solution[,2],
                            "x2" = IVP.solution[,3],
                            "x3" = IVP.solution[,4],
                            "x4" = IVP.solution[,5],
                            "x5" = IVP.solution[,6],
                            "x6" = IVP.solution[,7])
IVP.observation <- data.frame("time" = IVP.solution$time)
IVP.observation["y1"] <- NA
IVP.observation["y2"] <- NA
IVP.observation["y3"] <- NA

IVP.observation.unperturbed <- data.frame("time" = IVP.solution$time)
IVP.observation.unperturbed["y1"] <- NA
IVP.observation.unperturbed["y2"] <- NA
IVP.observation.unperturbed["y3"] <- NA

for (i in 1:length(time.sequence)) {
  IVP.observation[i,-1] <- system.PerturbedObservation(
    x = IVP.solution[i,-1],
    t = IVP.observation$time[i],
    params = system.parameters)
}
for (i in 1:length(time.sequence)) {
  IVP.observation.unperturbed[i,-1] <- system.ExpectedObservation(
    x = IVP.solution[i,-1],
    t = IVP.observation$time[i])
}
# <--- SOLVE IVP USING deSolve
# WRITE SAMPLED DATA INTO FILE --->

write( toJSON(IVP.observation), file = "C:/Users/kahl/Documents/Data/review_observation_json.txt")
write( toJSON(IVP.solution)   , file = "C:/Users/kahl/Documents/Data/review_truestate_json.txt")

write.table(IVP.solution, 'state.data.txt', sep = "\t")
write.table(IVP.observation,'observation.data.txt', sep = "\t")
write.table(IVP.observation.unperturbed, 'unperturbedobservation.data.txt', sep = "\t")
# <--- WRITE SAMPLED DATA INTO FILE