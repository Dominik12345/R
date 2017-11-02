
#' Class to handle all the equations
#'
odeEq <- setClass(
  #name of Class
  "odeEquations",
  slots = c(
    modelStr = "character",
    origEq = "character",
    measureFunction = "character",
    costateEq = "character",
    jacobian = "matrix",
    costFunction = "character",
    hamiltonian = "character",
    dynamicElasticNet = "logical"
  ),
  
  prototype = list(
    modelStr = character(0),
    origEq = character(0),
    measureFunction =  character(0),
    costateEq = character(0),
    costFunction = character(0),
    jacobian = matrix(list(),nrow = 2, ncol = 2),
    hamiltonian = character(0),
    dynamicElasticNet = FALSE
  ),
  
  validity = function(object) {
    if(!(is.character(object@origEq) && is.character(object@costateEq) && is.character(object@hamiltonian) && is.matrix(object@jacobian))) {
      return("The Equations have to be characters and the Jacobian a matrix.")
    }
    return(TRUE)
  }
)

setGeneric(name="setCostateEq",
           def = function(theObject,costate)
           {
             standardGeneric("setCostateEq")
           }
)

setMethod(f = "setCostateEq",
          signature = "odeEquations",
          definition = function(theObject,costate)
          {
            theObject@costateEq <- costate
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setOrigEq",
           def = function(theObject,origEq)
           {
             standardGeneric("setOrigEq")
           }
)

setMethod(f = "setOrigEq",
          signature = "odeEquations",
          definition = function(theObject,origEq)
          {
            theObject@origEq <- origEq
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setJacobian",
           def = function(theObject,jacobian)
           {
             standardGeneric("setJacobian")
           }
)

setMethod(f = "setJacobian",
          signature = "odeEquations",
          definition = function(theObject,jacobian)
          {
            theObject@jacobian <- jacobian
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setHamiltonian",
           def = function(theObject,hamiltonian)
           {
             standardGeneric("setHamiltonian")
           }
)

setMethod(f = "setHamiltonian",
          signature = "odeEquations",
          definition = function(theObject,hamiltonian)
          {
            theObject@hamiltonian <- hamiltonian
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "calculateCostate",
           def = function(theObject) {
             standardGeneric("calculateCostate")
           }
)

setMethod(f = "calculateCostate",
          signature = "odeEquations",
          definition = function(theObject)
          {
            tempList <- symbolicDiff(theObject)
            theObject@costateEq <- tempList$costate
            theObject@jacobian <- tempList$jacobian
            theObject@origEq <- tempList$origEq
            theObject@hamiltonian <- tempList$Hamilton
            return(theObject)
          }
)

setGeneric(name = "createModelEqClass",
           def = function(theObject,theModel) {
             standardGeneric("createModelEqClass")
           }
)

setMethod(f = "createModelEqClass",
          signature = "odeEquations",
          definition = function(theObject,theModel)
          {
            theObject@modelStr <- deparse(theModel,width.cutoff = 500)
            theObject@origEq <- getEquations(theModel)
            return(theObject)
          }
          
)

setGeneric(name = "setCostFunc",
           def = function(theObject,costFunction){
             standardGeneric("setCostFunc")
           }
)

setMethod(f = "setCostFunc",
          signature = "odeEquations",
          definition = function(theObject,costFunction)
          {
            theObject@costFunction <- getEquations(costFunction)
            theObject@dynamicElasticNet <- FALSE
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "setMeassureFunc",
           def = function(theObject,meassureFunc) {
             standardGeneric("setMeassureFunc")
           }
)

setMethod(f = "setMeassureFunc",
          signature = "odeEquations",
          definition = function(theObject,meassureFunc)
          {
            theObject@measureFunction <- getEquations(meassureFunc)
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "isDynElaNet",
           def = function(theObject) {
             standardGeneric("isDynElaNet")
           }
)

setMethod(f = "isDynElaNet",
          signature = "odeEquations",
          definition = function(theObject)
          {
            theObject@dynamicElasticNet <- TRUE
            return(theObject)
          }
)