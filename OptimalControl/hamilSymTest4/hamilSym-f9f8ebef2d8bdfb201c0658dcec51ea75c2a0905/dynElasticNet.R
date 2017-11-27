dynElasticNet <- function(alphaStep,armijoBeta,x0,parameters,times,alpha1,alpha2,measData,STD,modelFunc,measFunc,odeObj,optW,origAUC,maxIteration) {
  

  #' loading different libraries and checking the inputs
  #' loading
  library("deSolve")
  library('pracma')
  source("classOdeEquation.R")
  source("costate.R")
  source("stateHiddenInput.R")
  
  startOptW <- optW
  
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

  # setting alpha1 = 0 for this approach
  alphaDynNet <- list(a1 = alpha1, a2 = alpha2) # list of the alpha_1 and alpha_2 values



  
  # interpolation of the standard deviations
  if(is.null(STD)) {
    # if the standard deviations are missing the weights are set to 1 for every point of measurement, no weights will be applied
    measureData <- as.matrix(measureData)
    Q <- matrix(data = rep(1,length(measureData)), ncol = ncol(measureData))
    ## normalization of the weights based on the meassured data as suggested by Dominik Kahl
    Q <- Q / (measureData)
    Q[is.infinite(Q)] = 0
    Q[is.na(Q)] = 0
    interpQ <- apply(X = Q, MARGIN = 2, FUN = function(t) approx(x=measureTimes, y=t, xout = times))
    interpQ = do.call(cbind, lapply(interpQ, FUN = function(t) cbind(t$y)) )
    Q <- interpQ
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

    
    iter = 50
    alpha = alphaS
    arrayJ = rep(0,iter)
    
    for (i in 1:iter) {
      newW = oldW + alpha*gradStep
      #optW = as.nummeric(colSums(abs(newW)) > 0)
      
      
      input$optW <- optW
      input$w <- apply(X = newW, MARGIN = 2, FUN = function(x) approxfun(x = Tp, y = x, method = 'linear', rule=2))
      time <- seq(from = tInt[1], to = tInt[2], length.out = 300)
      solX <- ode(y = x0, times = time,func = hiddenInputState, parms = parameters, input=input)
      
      Tx <- solX[,1]
      x <- solX[,-1]

      yHat <- getMeassures(solX,measFunc)
      
      
      #interpolate the results
      input$interpX <- apply(X = x, MARGIN = 2, FUN = function(x) approxfun(x = Tx, y = x, rule=2, method = 'linear'))
      input$interpyHat <- apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))

      
      
      arrayJ[i] = costFunction(measureTimes,input,alphaDynNet)
      if ( i>1 && (arrayJ[i]>arrayJ[i-1]) && (arrayJ[i] < J[[currIter]])) {
        alpha = alphaS*stepBeta^(i-2)
        #print(paste0('i= ',i-1,' alpha= ',alpha,' J[W]=', arrayJ[i-1]))
        break
      }
      beta = stepBeta^(i)
      alpha = alphaS*beta
    }

    return(alpha)
  }
  
  
  #' see the estimated results for each iteration
  #' only used for debugging
  showEstimates <- function(measureTimes,AUCs,input, alpha2, J){
    y <- sapply(input$interpY, mapply, measureTimes)
    yhat <- sapply(input$interpyHat, mapply, measureTimes)
    
    minY = min(min(y),min(yhat))
    maxY = max(max(y),max(yhat))
    
    J <- unlist(J)
    J = J[J!=0]
    
    width = 2
    par(mfrow=c(1,3))
    barplot(unlist(AUCs[1,]), col = 'red', xlab = 'hidden inputs', main = 'AUC (a.u)')
    matplot(x = measureTimes, y = y, type = 'l', col = 4, xlab = '', ylab = '', ylim = c(minY,maxY), lwd = width)
    par(new=TRUE)
    matplot(x = measureTimes, y = yhat, type = 'l', col = 2, xlab = 't', ylab = 'y', ylim = c(minY,maxY), lwd = width)
    title(expression("y " * phantom("y") * " "), col.main = "blue")
    title(expression(phantom("y ") * hat('y')), col.main = "red")
    plot(J, type = 'l', xlab = 'iteration', ylab = 'J[w]', lwd = width)
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
  qInterp <- apply(X = Q, MARGIN = 2, FUN = function(x) approxfun(x =times, y = x, rule = 2, method = 'linear'))
  wInterp <- apply(X = w, MARGIN = 2, FUN = function(x) approxfun(x = Tx, y = x, rule = 2, method = 'linear'))
  
  # list of functions that approximate the data
  input <- list(optW= optW,interpX =xInterp,interpY = yInterp, interpyHat= yHatInterp, q=qInterp, w = wInterp)
  
  J <- list()
  J[[1]] <- costFunction(measureTimes,input,alphaDynNet)
  print(paste0('nominal cost J[w]= ',J[[1]]))
  
  offset <- 5
  usedAlphas <- rep(0,offset)
  
  if(missing(maxIteration)) {
    maxIter <- 50
  } else {
    maxIter <- maxIteration
  }
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
    #step = P - alpha2*w

    #' conj. gradient
    #' 
    if(i==1){
      gNeg = P - alpha2*w
      oldGrad = -gNeg
      step = gNeg
    }
    else {
      gNeg = P - alpha2*w
      newInt <- apply(X = -gNeg, MARGIN = 2, FUN = function(x) trapz(Tp, x^2))
      oldInt <- apply(X = oldGrad, MARGIN = 2, FUN = function(x) trapz(Tp, x^2))
      #betaStep =  newInt / oldInt
      betaStep = mldivide(oldInt,newInt)
      betaStep[is.nan(betaStep)] <- 0
      oldGrad = -gNeg

      #step = gNeg + t(t(step)*betaStep)
      step = gNeg + as.numeric(betaStep)*step

    }

    # calculating the stepsize
    if(i < offset) {
      alphaS = alphaStep
    } else {
      alphaS = max(usedAlphas)
      if(length(unique(usedAlphas))==1){
        alphaS = 4*alphaS
      }
    }
    stepBeta = armijoBeta
    alpha = getAlphaBacktracking(oldW,w,Q,y,step,J,i,alphaDynNet,alphaS,stepBeta,optW,parameters,tInt,Tp,measFunc,input,measureTimes)
    usedAlphas[mod(i,offset)+1] = alpha
    
    # calculate the new hidden inputs
    w = oldW + alpha*step
    
    
    #' CALCULATION OF THE TRAJEKTORIES THAT RESULTS FROM THE NEW HIDDEN INPUT
    inputState <- list()
    inputState$optW <- optW
    inputState$wInterp <- apply(X = w, MARGIN = 2, FUN = function(x) approxfun(x = Tp, y = x, method = 'linear', rule=2))
    
    

    solX <- ode(y = x0, times = times,func = hiddenInputState, parms = parameters, input=inputState)

    Tx <- solX[,1]
    x <- solX[,-1]
    

    yHat <- getMeassures(solX,measFunc)
    # interp
    input$interpX <- apply(X = x, MARGIN = 2, FUN = function(x) approxfun(x = Tx, y = x, rule=2, method = 'linear'))
    input$interpyHat <- apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
    input$w <- inputState$wInterp

    # #calculate the new cosT
    J[[i+1]] = costFunction(measureTimes,input,alphaDynNet)
    
    #AUCs <- sapply(X = input$w, FUN = function(x) trapzfun(f = x, a = t0, b = tf))
    tAUC <- seq(from = t0, to= tf, length.out =N)
    absW <- abs(sapply(input$w, mapply, tAUC))
    interpAbsW <- apply(X = absW, MARGIN = 2, FUN = function(x) approxfun(x = tAUC, y = x, rule=2, method = 'linear'))
    
    AUCs <- sapply(X = interpAbsW, FUN = function(x) trapzfun(f = x, a = t0, b = tf))
    #' output that shows the change of the costfunction and the aucs of the hidden inputs
    #print(paste0('change in costfunction: ',J[[i]]-J[[i+1]]))
    print(paste0('Iteration ',i,' J[w]=',round(J[[i+1]],8),'     change J[w]: ',round((1-abs(J[[i+1]]/J[[i]]))*100,2),' %  alpha=',alpha ))

    
    showEstimates(measureTimes,AUCs,input,alpha2,J)
    
    #' if the change in the cost function is smaller that epsilon the algorithmus stops
    #if (i == maxIter || abs(J[[i+1]]-J[[i]])< epsilon) {
    if (i == maxIter || abs(J[[i+1]]/J[[i]]) > 0.9975) {
      break
    }
  }
  
  
  if (missing(origAUC)) {
    origAUC <- AUCs
  }
  
  #' calculate the vector optW for the greedy approach
  #' the function calculated the combination of hidden inputs that could explain the deviation of the data
  #' the calculation is based on the aucs of the first optimisation
  
  greedySelection <- function(AUC, optW,origAUC) {
        orderAUCs <- order(-do.call(cbind,as.list(origAUC[1,])))
        tempOptW = rep(0,ncol(AUC))
      
        if(sum(unlist(origAUC[1,]))==sum(unlist(AUC[1,]))) {
          # if no aucs are given (first optimisation is running) -> select bigges AUCs
          tempOptW[orderAUCs[1]] = 1
        }
        else {
          tempOptW[orderAUCs[1:(sum(optW)+1)]] = 1
        }
        return(tempOptW)
  }
  
  #' calculate the error of the estimated measurements
  #' 
  rmse <- function(measureTimes,input){
    y <- sapply(input$interpY, mapply, measureTimes)
    yhat <- sapply(input$interpyHat, mapply, measureTimes)
    

    # feature scaling of the data
    y = apply(y, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
    yhat = apply(yhat, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
    
    return((colSums((yhat - y)^2))/nrow(yhat))
  }
  
  #' results are return in a list, because of the different data types
  
  results <- list()
  
  results$w <- w
  results$AUC <- do.call(cbind,AUCs[1,])
  results$optW <- greedySelection(AUCs, optW, origAUC)
  results$x <- solX
  results$y <- yHat
  results$rmse <- rmse(measureTimes,input)
  lastJ <- unlist(J)
  lastJ = lastJ[lastJ>0]
  lastJ = lastJ[length(lastJ)]
  results$J <- lastJ
  results$totalJ <- J
  
  
  return(results)
  
}