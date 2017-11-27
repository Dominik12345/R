greedyApproach <- function(alphaStep,Beta,alpha1, alpha2, x0, optW, times, measFunc, measData, std,
                           parameters, modelFunc, greedyBoolean) {
  graphics.off()
  
  source("classOdeEquation.R")
  source('symbolicDiff.R')
  source("createFunctions.R")
  source("dynElasticNet.R")
  
  #create the needed files
  odeEq <- odeEq()
  odeEq <- createModelEqClass(odeEq,modelFunc)
  odeEq <- setMeassureFunc(odeEq,measFunc)
  odeEq <- isDynElaNet(odeEq)
  odeEq <- calculateCostate(odeEq)
  createFunctions(odeEq)
  
  
  #' optW if knots are specific left out of optimisation, get the maximal estimated inputs
  iter <- (sum(optW))
  
  
  
  #' cross validation alpha
  #' tries to calculate the best value for alpha2 based on the best fit, without setting the value to low
  crossRes <- list()
  if(is.null(alpha2)) {
    steps <- 5
    error <- matrix(rep(0,2),ncol=2)
    colnames(error) <- c('alpha','MSE')
    for (i in 1:steps) {
      alpha1 = 0
      alpha2 = 1*10^(1-i)
      print('')
      print(paste0('Alpha 1=',alpha1,' Alpha2=',alpha2))
      crossRes[[i]] <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta,x0 = x0, optW = optW,
                               times=times, measFunc= testMessure, measData = y, STD = std,
                               alpha1 = alpha1, alpha2 = alpha2,
                               parameters = parameters, modelFunc = testModel, odeObj = odeEq,maxIteration=100)
      if (i==1){
        error[i,] = c(alpha2,mean(crossRes[[i]]$rmse))
      } else {
        error = rbind(error,c(alpha2,mean(crossRes[[i]]$rmse)))
      }

      print(error)
    }
    alpha2 = 1*10^(1-which.min(error[,2]))
    results <- crossRes[[which.min(error[,2])]] #use the estimated results of the crosssvalidation for saving time
  } else {
    #' Get initial values for the aucs
    #' start first esitmation
    results <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta, x0 = x0, optW = optW, 
                             times=times, measFunc= testMessure, measData = y, STD = std,
                             alpha1 = alpha1, alpha2 = alpha2,
                             parameters = parameters, modelFunc = testModel, odeObj = odeEq)
  }
  
  if (!greedyBoolean) {
    return(results)
  } 
  else {
    orgOptW <- optW <- results$optW
    orgAUC <- results$AUC
    
    barplot(orgAUC)
    print(results$rmse)
    
    optWs <- list()
    
    resAlg <- list()
    
    costError <- cbind(rep(0,length(optW)),rep(0,length(optW)))
    colnames(costError) <- c('sum(MSE)','cost')
    
    alphaStep = alphaStep*10
    
    
    for(i in 1:iter) {
      print('-----------------------------------------')
      print('selection done: starting new optimization')
      print('selected inputs:')
      print(which(optW > 0))
      optWs[[i]] <- optW
      resAlg[[i]] <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta, alpha1 = alpha1, alpha2 = alpha2,x0 = x0, optW = optW, 
                                   times=times, measFunc= testMessure, measData = y, STD = std, 
                                   parameters = parameters, modelFunc = testModel, odeObj = odeEq, origAUC = orgAUC)
      
      print(resAlg[[i]]$optW)
      
      
      costError[i,] = c(sum(resAlg[[i]]$rmse),resAlg[[i]]$J)
      
      
      
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
  
}
