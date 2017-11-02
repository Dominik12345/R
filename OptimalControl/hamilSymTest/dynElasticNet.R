dynElasticNet <- function(x0,parameters,times,alpha1,alpha2,measData,STD,modelFunc,measFunc,odeObj,optW,origAUC) {
  
  #' loading different libraries and checking the inputs
  #' loading
  library("deSolve")
  library('pracma')
  source("classOdeEquation.R")
  source("costate.R")
  
  # constructor for the neccessary functions
  checkInput <- function(odeObj) {
    if (class(odeObj)[1]!="odeEquations") {
      print("Please provide the neccessary functions for the dynamic elastic net approach.")
    }
    else {
      odeEq <- odeEq() # generate odeEquation S4-Object
      
    }
      
  }
  
  # Check for the optW vector, only states with an index of 1 will be calclated to enable the 'greedy' approach
  if (missing(optW)) {
    optW <- rep(1,length(x0))
  }
  
  #' Initial Options
  optInit <- sum(optW)  # Anzahl der zu optimierenden Inputs
  epsilon <- 0.1      # Abbruchkriterium der Armijo-Schrittweitenbestimmung

  N <- 100
  
  t0 <- times[1]
  tf <- tail(times, n=1)
  times <- seq(from = t0, to = tf, length.out = N)
  tInt <- c(t0,tf)
  
  
  measureTimes <- measData[,1] # select the time points of the given data
  measureData <- measData[,-1] # select the data
  if(missing(alpha2)) {
    alpha2 <- 0.2 
  }
  # setting alpha1 = 0 for this approach
  alpha1 <- 0
  alphaDynNet <- list(a1 = alpha1, a2 = alpha2) # list of the alpha_1 and alpha_2 values



  
  # interpolation of the standard deviations
  if(missing(STD)) {
    # if the standard deviations are missing the weights are set to 1 for every point of measurement, no weights will be applied
    Q <- matrix(data = rep(1,ncol(measureData)*length(times)), ncol = ncol(measureData))
  }
  else {
    # interpolating the given standard deviations to be used as weights in the costete equation
    interpSTD <- apply(X = STD, MARGIN = 2, FUN = function(t) approx(x = measData[,1], y = t, xout = times))
    interpSTD = do.call(cbind, lapply(interpSTD, FUN = function(t) cbind(t$y)))
    Q <- apply(X = interpSTD, MARGIN = 2, FUN = function(t)  (1/t^2)/length(t) )
  }

  #' get the trajectories of the nominal mondel without inputs
  #' using deSolve 'ode' function
  #' the first evaluation of the costfunction is based on the nominal model with hidden inputs constant at zero for all times
  if(all(is.null(names(x0)))) {
    names(x0) <- paste0(rep("x",length(x0)),1:length(x0))
  }
  solNominal <- as.data.frame(ode(y = x0, times = times, func = modelFunc, parms = parameters))
  
  # Tx the timepoints of the solution of the nominal model
  Tx <- solNominal[,1]
  # X trajectories of the states
  x <- solNominal[,-1]
  
  #initialize the hidden input vector
  w <- matrix(rep(0,nrow(x)*ncol(x)), nrow = nrow(x)) # is initialized as constant zero
  colnames(w) <- paste0(rep("w",ncol(w)),1:ncol(w))

  
  #' MEASURE OF THE CALCULATED TRAJEKTORIES
  getMeassures <- function(x,measFunc) {
    if(missing(measFunc)) {
      # if there is no function given for measurements it is assumed that all states are measurable 
      print('No meassurement function defined. Assuming all states are observable.')
      y <- x[,-1]
    } else {
      # combine the listed measurements to a matrix using do.call and cbind
      y <- do.call(cbind,measFunc(x[,-1]))
    }
    y = as.data.frame(cbind(x[,1],y))
    names(y)[-1] <- paste0(rep("y",ncol(y)-1),1:(ncol(y)-1)) 
    return(y)
  }
  
  # function to calculate the stepsize for the gradient descent
  getAlphaBacktracking <- function(oldW,W,Q,y,gradStep,J,currIter,alphaDynNet,alphaS,stepBeta,optW,para,tInt,Tp,measFunc,input,measureTimes) {
    iter = 20
    alpha = alphaS
    arrayJ = rep(0,iter)
    
    for (i in 1:iter) {
      newW = oldW + alpha*gradStep
      #optW = as.nummeric(colSums(abs(newW)) > 0)
      
      w <- apply(X = newW, MARGIN = 2, FUN = function(x) approxfun(x = Tp, y = x, method = 'linear', rule=2))
      time <- seq(from = tInt[1], to = tInt[2], length.out = 300)
      solX <- ode(y = x0, times = time,func = hiddenInputState, parms = parameters, input=w)
      
      Tx <- solX[,1]
      x <- solX[,-1]

      yHat <- getMeassures(solX,measFunc)
      
      
      #interpolate the results
      input$interpX <- apply(X = x, MARGIN = 2, FUN = function(x) approxfun(x = Tx, y = x, rule=2, method = 'linear'))
      input$interpyHat <- apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
      input$w <- w
      
      
      arrayJ[i] = costFunction(measureTimes,input,alphaDynNet)
      #print(paste0('i= ',i,' alpha= ',alpha,' J[W]=', arrayJ[i]))
      if ( i>1 && (arrayJ[i]>arrayJ[i-1]) && (arrayJ[i] < J[[currIter]])) {
        alpha = alphaS*stepBeta^(i-2)
        print(paste0('i= ',i-1,' alpha= ',alpha,' J[W]=', arrayJ[i-1]))
        break
      }
      beta = stepBeta^(i)
      alpha = alphaS*beta
    }
    return(alpha)
  }
  
  # cost function that is to be optimized
  costFunction <- function(measureTimes,input,alphaDynNet) {
    y <- sapply(input$interpY, mapply, measureTimes)
    yhat <- sapply(input$interpyHat, mapply, measureTimes)
    q <- sapply(input$q, mapply, measureTimes)
    w <- sapply(input$w, mapply, measureTimes)
    
    yCost <- list()
    # cost of the deviation of the calculated measurements to the given data
    for (i in 1:ncol(yhat)) {
      yCost$Start[[i]] = sum((yhat[1,i]- y[1,i]) * q[1,i] * (yhat[1,i]- y[1,i]))
      yCost$Middle[[i]] = sum((yhat[,i]- y[,i]) * q[,i] * (yhat[,i]- y[,i]))
      yCost$End[[i]] = sum((yhat[nrow(yhat),i]- y[nrow(y),i]) * q[nrow(q),i] * (yhat[nrow(yhat),i]- y[nrow(y),i]))
    }
    
    # cost for the inputs
    wCost <- list(L1 = 0, L2 = 0)
    for (i in 1:ncol(w)) {
      wCost$L1 = wCost$L1 + sum(abs(w[,i]))
      wCost$L2 = wCost$L2 + sum(abs(w[,i]^2))
    }
    
    #combining the costs
    cost = sum(yCost$Start)  + sum(yCost$Middle) + sum(yCost$End) + alphaDynNet$a1*wCost$L1 + alphaDynNet$a2*wCost$L2
    return(cost)
  }
  
  # calculate the measurements
  yHat <- getMeassures(solNominal,measFunc)
  
  # interpolation
  # linear approximation of the calculated values of x,y and yhat
  xInterp <- apply(X = x, MARGIN = 2, FUN = function(x) approxfun(x = Tx, y = x, rule=2, method = 'linear'))
  yInterp <- apply(X = measData[,-1], MARGIN = 2, FUN = function(x) approxfun(x = measData[,1], y = x, rule=2, method = 'linear'))
  yHatInterp <- apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
  qInterp <- apply(X = Q, MARGIN = 2, FUN = function(x) approxfun(x = times, y = x, rule = 2, method = 'linear'))
  wInterp <- apply(X = w, MARGIN = 2, FUN = function(x) approxfun(x = Tx, y = x, rule = 2, method = 'linear'))
  
  # list of functions that approximate the data
  input <- list(interpX =xInterp,interpY = yInterp, interpyHat= yHatInterp, q=qInterp, w = wInterp)
  
  J <- list()
  J[[1]] <- costFunction(measureTimes,input,alphaDynNet)
  print(J)


  
  maxIter <- 100
  for (i in 1:maxIter) {
    
    #' CALCULATION OF THE HIDDEN INPUTS BASED ON THE COSTATES AND OLD INPUTS
    #' #####################################################################
    #' calculation of the costates
    costateStart <- rep(0,ncol(x))
    timesCostate <- seq(from = tf, to= t0, length.out =N)
    
    solCostate <- ode(y = costateStart, times = timesCostate, func = costate, parms = parameters, input=input)
    solCostate = solCostate[nrow(solCostate):1,]

    Tp = solCostate[,1]
    P = solCostate[,-1]
    
    
    
    oldW = w
    # gradient step
    step = P - alpha2*w
    
    # calculating the stepsize
    alphaS = 0.01
    stepBeta = 0.5
    alpha = getAlphaBacktracking(oldW,w,Q,y,step,J,i,alphaDynNet,alphaS,stepBeta,optW,parameters,tInt,Tp,measFunc,input,measureTimes)
    
    # calculate the new hidden inputs
    w = oldW + alpha*step
    
    
    
    #' CALCULATION OF THE TRAJEKTORIES THAT RESULTS FROM THE NEW HIDDEN INPUT
    wInterp <- apply(X = w, MARGIN = 2, FUN = function(x) approxfun(x = Tp, y = x, method = 'linear', rule=2))
    solX <- ode(y = x0, times = times,func = hiddenInputState, parms = parameters, input=wInterp)
    Tx <- solX[,1]
    x <- solX[,-1]

    yHat <- getMeassures(solX,measFunc)
    # interp
    input$interpX <- apply(X = x, MARGIN = 2, FUN = function(x) approxfun(x = Tx, y = x, rule=2, method = 'linear'))
    input$interpyHat <- apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
    input$w <- wInterp

    # #calculate the new cosT
    J[[i+1]] = costFunction(measureTimes,input,alphaDynNet)
    
    AUCs <- sapply(X = input$w, FUN = function(x) trapzfun(f = x, a = t0, b = tf))

    #' output that shows the change of the costfunction and the aucs of the hidden inputs
    print(paste0('change in costfunction: ',J[[i]]-J[[i+1]]))
    print(do.call(cbind,AUCs[1,]))
    
    #' if the change in the cost function is smaller that epsilon the algorithmus stops
    if (i == maxIter || abs(J[[i+1]]-J[[i]])< epsilon) {
      
      break
    }
  }
  
  
  #' calculate the vector optW for the greedy approach
  #' the function calculated the combination of hidden inputs that could explain the deviation of the data
  #' the calculation is based on the aucs of the first optimisation
  greedySelection <- function(AUC, optW,origAUC) {
    if (sum(optW) != length(optW)) {
      orderAUCs <- order(-do.call(cbind,origAUC[1,]))
      indSel <- orderAUCs[1:(sum(optW)+1)]
      optW = rep(0,ncol(AUC))
      optW[indSel] = 1
    } 
    else {
      orderAUCs <- order(-do.call(cbind,origAUC[1,]))
      indSel <- orderAUCs[1]
      optW = rep(0,ncol(AUC))
      optW[indSel] = 1
    }
    return(optW)
  }
  
  #' results are return in a list, because of the different data types
  
  results <- list()
  if(sum(optW) == length(optW)) {
    results$origAUC = AUCs
    print(results$origAUC)
  }
  
  results$w <- w
  results$AUC <- do.call(cbind,AUCs[1,])
  results$optW <- greedySelection(AUCs, optW, AUCs)
  results$x <- solX
  results$y <- yHat
  
  
  return(results)
  
}