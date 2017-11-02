#' AUX FUNCTIONS
#' 
#' a set of functions that are simply there for a better workflow
#' 

## Check on which machine is worked on and set wd corretly
setWd <-function() {
  system <- R.Version()$system
  if(grepl("linux",system)) {
    setwd("/home/newmi/hamilSym")
  }
  else {
    setwd('D:/R-Programmierung/symbolic')
  }
}

## Check if needed packages are installed. If packages are not present, R will install them
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

strMatrixMulti <- function(M1,M2){ 
  if(ncol(M1) != nrow(M2)) {
    print('multiplication not possible. Dimensions dont match')
  }
  else {
    colM1 <- ncol(M1)
    colM2 <- ncol(M2)
    rowM1 <- nrow(M1)
    rowM2 <- nrow(M2)
    
    ret <- matrix(rep(0,rowM1*colM2), ncol=colM2)
    m <- rowM1
    n <- colM2
    for (r in 1:m) {
      for (c in 1:n) {
        multi <- paste(M1[r,],M2[,c],sep = "*")
        ret[r,c] = paste(multi,collapse = "+")
      }
    }
  }
  return(ret)
}
