% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seeds.R
\docType{package}
\name{seeds-package}
\alias{seeds}
\alias{seeds-package}
\title{seeds: Estimate hidden inputs using the method of the dynamic elastic net}
\description{
A package combining two different algorithms for calculating hidden
inputs. These hidden inputs can be interpreted as a controll that tries to minimize the
discrepancies between a given model and a given measurements. The method is
also called the 'dynamic elastic net'.
}
\details{
Details

The algorithm calcultes the needed equations using the \code{\link{Deriv}}
function of the \pkg{Deriv} package. The process is implemented through the use
of the S4 class \code{\link{odeEquations-class}}. 

The conjugate gradient based algorithm uses a greedy algorithm to estimate a 
sparse controll that tries to minimize the discrepancies between a given 
'nominal model given the measurements (e.g from an experiment). The algorithm
the \code{\link{ode}} uses \pkg{deSolve} to calculate the hidden inputs w 
based on the adjoint equations of the ODE-System. 

The adjoint equations are calculated using the \code{\link{ode}} function of the 
\pkg{deSolve} package. For the usage of the algorithm please look into the 
examples and documentation given for the functions 
\describe{
   \item{\code{\link{greedyApproach}}}{a greedy algorithm to calculate a sparse control}
   \item{\code{\link{dynElasticNet}}}{a conjugate gradient optimal control algorithm with L2 regulation} 
}
}
\references{
\strong{Benjamin Engelhardt, Holger Fröhlich, Maik Kschischo} 
Learning (from) the errors of a systems biology model, \emph{Nature Scientific Reports},
 6, 20772, 2016 \url{https://www.nature.com/articles/srep20772}
}
\author{
\strong{Maintainer}: Tobias Newmiwaka \email{tobias.newmiwaka@gmail.com}

}
