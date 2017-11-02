greedyApproach <- function(alphaStep, x0, optW, times, measFunc, measData,  parameters, modelFunc, odeObj, greedyBoolean) {
  graphics.off()
  
  #' Get initial values for the aucs
  #' start first esitmation
  source("dynElasticNet.R") # load the algorithm
  print(optW)
  results <- dynElasticNet(alphaStep = alphaStep,x0 = x0, optW = optW, times=times, measFunc= testMessure, measData = y,  parameters = parameters, modelFunc = testModel, odeObj = odeEq)
  
  orgOptW <- optW <- results$optW
  orgAUC <- results$AUC
  
  barplot(orgAUC)
  print(results$rmse)
  
  optWs <- list()
  
  resAlg <- list()
  
  costError <- cbind(rep(0,length(optW)),rep(0,length(optW)))
  colnames(costError) <- c('sum(RMSE)','cost')
  
  for(i in 1:(length(optW)-1)) {
    print('-----------------------------------------')
    print('selection done: starting new optimization')
    cat('selected hidden inputs: ', optW)
    optWs[[i]] <- optW
    resAlg[[i]] <- dynElasticNet(alphaStep = alphaStep,x0 = x0, optW = optW, times=times, measFunc= testMessure, measData = y,  parameters = parameters, modelFunc = testModel, odeObj = odeEq, origAUC = orgAUC)
    
    print(resAlg[[i]]$optW)

    
    costError[i,] = c(sum(resAlg[[i]]$rmse),resAlg[[i]]$J)
    print(costError)
    
    
    
    if(i > 1 && ( costError[i,1] > costError[i-1,1])  ) {
      print('hidden inputs on knots:')
      print(which(optWs[[i-1]] %in% 1))
      break
    }
    optW <- resAlg[[i]]$optW
  }
  
  minY = min(min(measData[,-1]),min(resAlg[[i-1]]$y[,-1]))
  maxY = max(max(measData[,-1]),max(resAlg[[i-1]]$y[,-1]))
  
  
  par(mfrow=c(2,2))
  barplot(resAlg[[i-1]]$AUC)
  matplot(x = measData[,1], y = measData[,-1], type = 'l', col = 4, xlab = '', ylab = '', ylim = c(minY,maxY))
  par(new=TRUE)
  matplot(x = resAlg[[i-1]]$y[,1], y = resAlg[[i-1]]$y[,-1], type = 'l', col = 2, xlab = 't', ylab = 'y', ylim = c(minY,maxY))
  matplot(x = resAlg[[i-1]]$y[,1], y = resAlg[[i-1]]$w, type = 'l', col = 2, xlab = 't', ylab = 'w')
  plot(x= 1:length(unlist(resAlg[[i-1]]$totalJ)),y= unlist(resAlg[[i-1]]$totalJ), type = 'l', xlab = 'Iteration', ylab = 'J[w]')
  
  
  
  return(resAlg[[i-1]])
}
