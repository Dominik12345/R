createFunctions <- function(odeEq){ 
  
  trim <- function (x)  sub("^\\s+", "", x) 
  trimSpace <- function (x) gsub("\\s", "",x)
  
  formatEqs <- function(modelEq){
    return(gsub("(d[x,y]*)\\[([0-9]*)\\]=","",modelEq))
  }

  createCostate <- function(modelODE){
    file.create('costate.R')
    fileCostate <- file('costate.R')
    writeLines(modelODE,fileCostate)
    close(fileCostate)
  }
  
  createStateHiddenInput <- function(modelODE){
    file.create('stateHiddenInput.R')
    fileState <- file('stateHiddenInput.R')
    writeLines(modelODE,fileState)
    close(fileState)
  }
  
  formatFuncString <- function(odeEq,funcType) {
    formatStr <- gsub(pattern = "\\[|[1-9]|\\]",replacement = "", strsplit(odeEq@origEq,"=")[[1]][1])
    if(funcType=="costate"){
      stringOde <- odeEq@costateEq
      
      #generate new function 'header'
      
      #get additional parameters 
      paraInput <- paste0(gsub("^\\s+|\\s+$|)", "",strsplit(odeEq@modelStr[1],split = ",")[[1]][-(1:2)]),collapse = ",") 
      #check if cost function needs additional inputs
      addInputs <- getAddInputs(stringOde)
      #check if additional variables are defined in model-function
      formatStr = sub(pattern = "d",replacement = "",formatStr)
      addInputs = trimAddInputs(odeEq,addInputs,paraInput, formatStr)

      funcStartStr <- paste0(funcType," <-function(t,p,parameters,input) {")
      # funcStartStr <- append(funcStartStr,paste0("\twith (as.list(",paste(paraInput, collapse = ","),") {\n"),after = length(funcStartStr)+1)
      funcStartStr <- append(funcStartStr,paste0("\twith (as.list(",paraInput,"), {\n"),after = length(funcStartStr)+1)
      
      # get additional operations, like setting parameters
      addiOperUntrimed <- odeEq@modelStr[!grepl("\\{|\\}",odeEq@modelStr)]
      addiOper <- append(addiOperUntrimed[!grepl("function|list|d[x,y]",addiOperUntrimed)],"\n")
      
      if(length(addiOper)!=0){
        funcStartStr = append(funcStartStr,addiOper,after = length(funcStartStr)+2)
      }
      
      funcStartStr = append(funcStartStr,"\t\toptW <- input$optW", after = length(funcStartStr)+1)
      
      if(sum(as.numeric(grepl("yhat",addInputs)))>0) {
        dynNetInterpolate <- c(    "\t\tx <- sapply(input$interpX, mapply, t)",
                                   "\t\ty <- sapply(input$interpY, mapply, t)",
                                   "\t\tyhat <- sapply(input$interpyHat, mapply, t)",
                                   "\t\tq <- sapply(input$q, mapply, t)\n")
        
        funcStartStr = append(funcStartStr,dynNetInterpolate, after = length(funcStartStr)+2)
      }

      
      
      #generate new function wrap up
      toList <- paste0(unlist(strsplit(x = stringOde, split = "="))[seq(1, 2*length(stringOde), by = 2)],collapse = ",")
      funcEndStr <- paste0("\n\t\t\t\tlist(c(",toList,"))\n","\n  })\n}")
      stringOde <- paste0("\t\t\t\t",stringOde)
      topMid <- append(funcStartStr,stringOde, after = length(funcStartStr)+2)
      res <- append(topMid,funcEndStr, after = length(topMid)+1)
    }
    else {
      stringOde <- odeEq@origEq

      funcStartStr <- paste(funcType," <- function(t, x, parameters, input) { \n with (as.list(parameters),{ \n")
      funcStartStr = append(funcStartStr,"\t\toptW <- input$optW", after = length(funcStartStr)+1)
      
      addiOperUntrimed <- odeEq@modelStr[!grepl("\\{|\\}",odeEq@modelStr)]
      addiOper <- append(addiOperUntrimed[!grepl("function|list|d[x,y]",addiOperUntrimed)],"\n")

      
      
      dynNetInterpolate <- "\t\tw <- sapply(input$w, mapply, t)\n"
      funcStartStr = append(funcStartStr,dynNetInterpolate, after = length(funcStartStr)+1)
      if(length(addiOper)!=0){
        funcStartStr = append(funcStartStr,addiOper,after = length(funcStartStr)+2)
      }
      stringOde <- paste0("\t\t",paste0(stringOde,paste0("+ optW[",1:length(stringOde),"] *w[",1:length(stringOde),"]")))
      funcMiddleStr = append(funcStartStr,stringOde, after = length(funcStartStr)+2)
      
      toList <- paste0(unlist(strsplit(x = stringOde, split = "="))[seq(1, 2*length(stringOde), by = 2)],collapse = ",")
      funcEndStr <- paste0("\n\t \tlist(c(",gsub("\t","",toList),"))","\n  })\n}")
      res <- append(funcMiddleStr,funcEndStr, after = length(funcMiddleStr)+1)
    }

    return(res)
  }
  
  getAddInputs <- function(stringVec) {
    arrayEle <- gsub('([a-z])\\[[0-9]\\]',stringVec, replacement = "\\1") # find elements with index
    arrayEleTrim <- gsub('[\\*, \\/]', arrayEle, replacement = " ")  #remove mathematical symbols
    arrayEleTrim = gsub('\\[|\\]|[+,-]|[=]|\\(|\\)',arrayEleTrim,replacement = " ") # get rid of the brackets
    arrayEleTrim = gsub(pattern = '[0-9]',x = arrayEleTrim, replacement = "")
    
    arrayUniqueEntries <- unique(unlist(strsplit(arrayEleTrim,split = " "))) #get unique entries
    arrayUniqueEntries <-arrayUniqueEntries[arrayUniqueEntries != ""] #remove empty entries
    arrayUniqueEntries = arrayUniqueEntries[!grepl("dp|p", arrayUniqueEntries)] # subsetting so no mathes for the costate parameter p are in
    
    return(arrayUniqueEntries)
  }
  
  # function to check for setting variables
  trimAddInputs <- function(odeEq,addInputs,paraInput, formatStr) {
    
    # get the definitions from the head of the eq
    linesDef <- tolower(trimSpace(trim(odeEq@modelStr[grepl("<-|[=]",odeEq@modelStr)]))) #models with
    eq <- trimSpace(odeEq@origEq)
    # subsetting to get rid of model equation
    defines <- linesDef[!linesDef %in% eq]
    # split the variable 'devines' to get the defined variables
    #definedPara
    defines = unlist(strsplit(defines,split = "<-|="))
    definedPara <- defines[seq(from = 1, to = 2*length(defines), by = 2)]
    paraInputEq <- defines[seq(from = 2, to = 2*length(defines), by = 2)]
    
    definedPara = definedPara[!is.na(definedPara)]
    paraInputEq = paraInputEq[!is.na(paraInputEq)]
    
    
    definedPara = unique(gsub(pattern = "[0-9]*",replacement = "", x = definedPara))
    paraInputEq = unique(gsub(pattern = "\\[[a-z]*[0-9]*\\]",replacement = "", x = paraInputEq))

    addInputs <- addInputs[!addInputs %in% definedPara]
    addInputs = append(paraInputEq,addInputs)

    if(sum(as.numeric(grepl(pattern = formatStr, addInputs)))>0) {
      
      return(addInputs)
    } 
    else {
      return(append(addInputs,formatStr))
    }
    #return(addInputs)
  }
  
  
  wrapperCreateFunc <- function(odeEq) {
    createCostate(formatFuncString(odeEq,"costate"))
    if(odeEq@dynamicElasticNet){
      createStateHiddenInput(formatFuncString(odeEq,"hiddenInputState"))
    }
  }

  wrapperCreateFunc(odeEq)
  
}