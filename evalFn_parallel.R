library(parallel) 
library(doParallel)
library(foreach)
library(iterators)
nCores <- 4  
registerDoParallel(nCores) 
# to addd
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

singleEval <- function(singleGene, X, y, type, criterion, criFun, family){
  if(type == "lm")
    fit <- lm(y~as.matrix(X[,which(singleGene != 0)]))
  if(type == "glm")
    fit <- glm(as.vector(y)~as.matrix(X[,which(singleGene != 0)]), family)

  if(is.null(criFun)){    # Dont have their own criterion function written
    criFunBuilt <- eval(parse(text = criterion))
    criValue <- criFunBuilt(fit)
  }
  else 
    criValue = try(criFun(fit), silent = TRUE)   # use the function inputted by the user; CHECK FOR ERRORS?
    if(!is.null(attributes(criValue)$class))
      if(attributes(criValue)$class == "try-error")
        stop(cat(paste("criFun is not compatible. The following error occured:\n", geterrmessage())))
    if(length(criValue)!=1)
      stop("Dimension of output for criFun greater than 1.")
    if(!is.numeric(criValue)!=1)
      stop("Output for criFun is not numeric.")
  return(criValue)
}
  
  
evalFunction2 <- function(currentGenePool, popSize, type = "lm", family = "gaussian", 
                          criterion = "AIC", criFun = NULL){
  if(criterion!="AIC" & criterion!= "BIC")  ### ADD MORE?
    stop(paste(criterion, "is not a valid criterion. Please use AIC or BIC."))
  
  if(!is.null(criFun) & !is.function(criFun))
    stop("criFun input is not a function.")
  
  if(type != "lm" & type != "glm")
    stop("Regression must be of type 'lm' or 'glm'")
  
  if(family == "binomial" & length(unique(na.omit(y)))!= 2)
    stop("Logistic regression requires 'y' to be binary")
  
  geneLength <- dim(currentGenePool)[2]
  result <- rep(NA, popSize)
   
  result <- foreach(i = 1:popSize, .combine = c)  %dopar% {
    criValue <- singleEval(currentGenePool[i,], X, y, type, criterion, criFun, family)
    return(criValue)
  }

  obj <- rbind(result, rank(result), rank(-result)/sum(1:popSize))
  row.names(obj) <- c(criterion, "ranks", "samplingProbs")
  return(obj)
}

#evalFunction2(genePool, popSize = 4, type = "lm", criterion = "AIC", criFun = NULL)

