library(deSolve)

#get the system equations and initial value
source('SystemEquation.R')

# SOLVE IVP USING deSolve --->

#   must be defined in SystemEquation:
#       state.initial 
#       parameters.updatestate 

IVP.Function <- function(t, state, parameters) {
  with(as.list(c(state,parameters)),{
    dx <- system.Dynamic(x = state, t = t, params = parameters ,delta.t = time.delta)
    return(list(dx) )
  })
}
sample.state <- c("x1" = system.parameters[["x1.0"]], "x2" = system.parameters[["x2.0"]])
IVP.solution <- ode(y = sample.state, 
                    times = time.sequence, 
                    func = IVP.Function, 
                    parms = system.parameters)

# convert into dataframe
IVP.solution <- data.frame( "time" = IVP.solution[,1], "x1" = IVP.solution[,2], 
                            "x2" = IVP.solution[,3])
# <--- SOLVE IVP USING deSolve

# SAMPLE OBSERVATION DATA --->
observation.data <- data.frame("time" = IVP.solution$time)
observation.data["y1"] <- NA
observation.data["y2"] <- NA
for (i in 1:length(time.sequence)) {
  observation.data[i,-1] <- system.Observation(x = IVP.solution[i,-1],
                                               t =  observation.data$time[i],
                                               params = system.parameters)
}
# <--- SAMPLE OBSERVATION DATA

# WRITE DATA INTO TXT --->
write.table(IVP.solution ,'state.data.txt', sep="\t")
write.table(observation.data, 'observation.data.txt',sep = "\t")
# <--- WRITE DATE INTO TXT '
