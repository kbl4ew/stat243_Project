##### Evaluation Function #####
##### Descritption ######
### Here we create two evaluation functions - evalLm and evalGlm
### that allows us evaluate the fitness of the current generation. ###

##### Input #####
### genePool: a matrix holding the current individuals of size population size by 
### X: the independent variables (training data)
### y: the dependent varaibles (training data)
### family: indicating the model for the general linear model function
### criterion: a character string specifying the function used as the evaluation criteria(default: "AIC")
### criFun: the user-provided function to use as evaluation criteria if "criterion" is not AIC or BIC

##### Output #####
### result: a matrix which holds the value from the criterion function & rank (assume
###         that the smaller criterion value, the smaller the rank, the higher weights)
### First row of result: evaluation criterion value
### Second row of result rank
#library(data.table)

evalLm <- function(genePool, X, y, criterion = "AIC",criFun = NA){
  
  popSize <- dim(genePool)[1] 
  geneLength <- dim(genePool)[2]
  result <- matrix(nrow = 1, ncol = popSize)
  rownames(result) <- c("criterion value")
  ### a for loop to iterate across each individual to calculate 
  for (i in 1:popSize){
    fit <- lm(y~as.matrix(X[,which(genePool[i,] != 0, arr.ind = T)]))
    if (criterion == "AIC"){
      criValue = AIC(fit)
    } else if (criterion == "BIC") {
      criValue = BIC(fit)
    } else {
      criValue = criFun(fit)
    }
    result[1,i] = criValue
  }
#   DT <- as.data.table(t(result))
#   DT[, AICRank := rank(AIC, ties.method = "first")]
#   return(t(DT))
  
  obj <- rbind(result, rank(result), rank(result)/sum(1:popSize))
  row.names(obj) <- c(criterion, "ranks", "samplingProbs")
  return(obj)
}


#evalLm()
#dim(mtcars)

#evalLm(test, X, y)
#rank(AICVals) / sum(1:length(AICVals)

evalGlm <- function(genePool, X, y,family = "gaussian", criterion = "AIC", criFun){
  popSize <- dim(genePool)[1] 
  geneLength <- dim(genePool)[2]
  result <- matrix(nrow = 1, ncol = popSize)
  rownames(result) <- c("criterion value")
  
  for (i in 1:popSize){
    #print("entered loop")
    fit <- glm(as.vector(y)~as.matrix(X[,which(genePool[i,] != 0, arr.ind = T)]),family)
    #print("got through this iteration")
    if (criterion == "AIC"){
      criValue = AIC(fit)
    } else if (criterion == "BIC") {
      criValue = BIC(fit)
    } else {
      criValue = criFun(fit)
    }
    result[1,i] = criValue
  }
  
  #   DT <- as.data.table(t(result))
  #   DT[, AICRank := rank(AIC, ties.method = "first")]
  #   return(t(DT))
  
  obj <- rbind(result, rank(result), rank(result)/sum(1:popSize))
  row.names(obj) <- c(criterion, "ranks", "samplingProbs")
  return(obj)
}

evalFunction <- function(currentGenePool, type, criterion, family = NA, criFun = NULL){
  if(type == "lm"){
    #print('lm flow')
    return(evalLm(currentGenePool, X, y, criterion = "AIC",criFun))
  }
  else if (type == "glm"){
    #print('glm flow')
    return(evalGlm(currentGenePool, X, y, family, criterion, criFun))
  }
  ### TO be fixed here
  else{
    return(c())
  }
  
}










