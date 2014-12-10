############################ Project: Main Code #################################
##### Description: #####
### This is the main code of the generic algorithm that attempts to solve the variable selection problem.

##### Input: #####
### Training Set: (X, y) (User Defined) #####
### genePool: size/# of individuals of the entire population
### mutation_prob: probability of mutation within each iteration
### crossover_prob: probability of crossover within each iteration

##### Output: #####
### result <- the best model based on the evaluation criterion
#X <- mtcars[,2:11]
#y <- mtcars[,1]



############################ Project: Main Code #################################
##### Description: #####
### This is the main code of the generic algorithm that attempts to solve the variable selection problem.

##### Input: #####
### Training Set: (X, y) (User Defined) #####
### genePool: size/# of individuals of the entire population
### mutation_prob: probability of mutation within each iteration
### crossover_prob: probability of crossover within each iteration
### max_iterations: maximum iterations to go through
### min_iterations

##### Output: #####
### result <- the best model based on the evaluation criterion
### Test Case (To be deleted)

##### Auxilary Functions #####
##### popInitialize #####
##### Input #####
# 1. popSize: The population size in each generation
# 2. geneLength: the number of genes in the chromosome
# 3. zeroToOneRatio: the change from a zero to one for mutations and initialization. (This option allows us to control the number of set genes in the chromosome. Since we are conducting variable selection, this parameter may be higher than usual.)

##### Output #####
# The output is a matrix of size (popSize, geneLength) with the values initiated
# We do not want to any individual with all zeros hence we guarantee that at least one element is 1.

##### Reference: genalg package
# The code is an adoption of the basic generic algorith implementation in the genalg package in R

##### Implementation #####
popInitialize <- function(popSize = 50, geneLength = 0, zeroToOneRatio = 1){
  
  if(geneLength == 0){
    stop("geneLength cannot be zero.")
  }

  pop <- matrix(nrow = popSize, ncol = geneLength);
  
  ##### Randomly initialize the first generation #####
  for (child in 1:popSize){
    pop[child, ] = sample(c(rep(0, zeroToOneRatio), 1), geneLength, replace = TRUE);
    
    while(sum(pop[child,]) == 0){
      pop[child, ] = sample(c(rep(0, zeroToOneRatio), 1), geneLength, replace = TRUE);   
      
    }
  }
  return(pop)
}

##### Auxilary Function 2: EvalFunction #####
#library(parallel) 
#library(doParallel)
#library(foreach)
#library(iterators)
#nCores <- 4  
#registerDoParallel(nCores) 
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
    fit <- lm(y~X[,which(singleGene != 0)])
  if(type == "glm")
    fit <- glm(as.vector(y)~X[,which(singleGene != 0)], family)
  
  if(is.null(criFun)){    # Dont have their own criterion function written
    criFunBuilt <- eval(parse(text = criterion))
    criValue <- criFunBuilt(fit)
  }
  else {
     criValue = try(criFun(fit), silent = TRUE)   # use the function inputted by the user; CHECK FOR ERRORS?
     if(!is.null(attributes(criValue)$class))
       if(attributes(criValue)$class == "try-error")
         stop(cat(paste("criFun is not compatible. The following error occured:\n", geterrmessage())))
     if(length(criValue)!=1)
       stop("Dimension of output for criFun greater than 1.")
     if(!is.numeric(criValue)!=1)
       stop("Output for criFun is not numeric.")
     return(criValue)
    return(0)
  }
  return(criValue)
}


evalFunction <- function(X, y, currentGenePool, popSize, type = "lm", family = "gaussian", criterion = "AIC", criFun = NULL){
  if((popSize %% 2)!= 0){
    message("Warning: The size of the population has been rounded to the largest even numer")
    popSize <- popSize + 1;
  }
  
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

##### Auxilary Function: UpdateSam #####
updateSamp <- function(x, popSize, weights){
  pairs = matrix(NA, ncol = popSize/2, nrow = 2) # form sampled pairs
  
  set.seed(0) # take out later
  for(i in 1:popSize/2){
    pairs[,i] = sample(1:popSize, 2, prob = weights)
  }
  pairs # so 2 and 6 are paired up; 4 and 5 paired up; 6,5 paired
  xSamp = x[as.vector(pairs), ]
  return(xSamp)
}

##### Auxilary Function: Crossover #####
crossover <- function(v1, v2, geneLength, crossRate){
  crossBool = sample(c(TRUE, FALSE), 1, prob = c(crossRate, 1-crossRate))
  
  if(crossBool){
    cut = sample(geneLength-1, 1)
    new1 = c(v1[1:cut], v2[(cut+1):geneLength])
    new2 = c(v2[1:cut], v1[(cut+1):geneLength])
    
    if(sum(new1) == 0 | sum(new2) == 0) # if either ended up with only zeros
      return(rbind(v1,v2))
    else
      return(rbind(new1, new2))
  }
  else
    return(rbind(v1,v2)) # return them unchanged
}

##### Auxilary Function: Mutation #####
mutation <- function(v1, v2, mRate){
  mLoci = which(v1==v2)
  len = length(mLoci)
  
  # T/F: mutate or not
  mBool1 = sample(c(TRUE, FALSE), len, replace = TRUE, prob = c(mRate, 1-mRate))
  mBool2 = sample(c(TRUE, FALSE), len, replace = TRUE, prob = c(mRate, 1-mRate))
  
  if(sum(mBool1) == 0 & sum(mBool2) == 0) 
    return(rbind(v1,v2)) # return v1 v2 and dont mutate
  
  else{
    v1Copy = v1
    v2Copy = v2
    v1Copy[mLoci][mBool1] <- as.numeric(!v1Copy[mLoci][mBool1])
    v2Copy[mLoci][mBool2] <- as.numeric(!v2Copy[mLoci][mBool2])
    
    if(sum(v1Copy) == 0) v1Copy <- v1
    if(sum(v2Copy)== 0) v2Copy <- v2
    
    return(rbind(v1Copy,v2Copy))    
  }
}


##### Clear out later #####
#X <- mtcars[,2:11]
#y <- mtcars[,1]

select <- function(X = NULL, y = NULL, popSize = 200, criterion = "AIC", type = "lm", family = "gaussian", criFun = NULL, max_iterations = 500, min_iterations = 50, crossRate = 0.95, mRate = 0.001, zeroToOneRatio = 1){
  ##### Defense coding #####
  #X <- as.matrix(X);
  #y <- as.vector(y);
  if((popSize%%2)!=0){
    #warning("The number of models has ")
    #print("Warning: The number of models has been incremented to the nearest even number")
    warning("The number of models has been incremented to the nearest even number")
    popSize <- popSize + 1
  }

  
  if(is.null(X)){
    stop("Please provide the predictors! Exiting from the function")
  }
  
   if(is.null(y)){
     stop("Please provide the independent variable/outcome! Exiting from the function")
  }
  library(parallel)
  library(doParallel)
  library(foreach)
  library(iterators)
  nCores <- 4  
#   checkLoaded <- try(registerDoParallel(nCores), silent = TRUE)
#   if(as.character(attributes(checkLoaded)[1] == "try-error")){
#     stop(cat("Please load the following packages before running select:\n1) parallel\n2) doParallel\n3) foreach\n4) iterators")) 
#   }
  ##### Beginning of the generic algorithm #####
  geneLength <- dim(X)[2];
  ##### Initializing the first generation of individuals/ models
  initialPopulation <- popInitialize(popSize, geneLength, zeroToOneRatio);
  currentGenePool <- initialPopulation;
  ### Calculating the sampling probabilities for the first generations of individuals/models
  #samplingProb <- evalFunction(currenGenePool, type, criterion, family, criFun)[3,]
  samplingProb <- evalFunction(X, y, currentGenePool, popSize, type, family, criterion, criFun)[3,];
  avgCriterion <- mean(evalFunction(X, y, currentGenePool, popSize, type, family, criterion, criFun)[1,]);
  ### While loop to handle convergence/ exceedance of min iteration/ capped by max iteration
  #iter = 0;
  ### Condition to be satisfied 
  ### if iter < min_iteration
  #while((iter <= min_iterations)&& !(iter >= max_iterations))
  for(i in 1:max_iterations){
    # really we will have predetermined # of iterations
    #xSamp <- updateSamp(x, popSize, weights)
    geneSample <- updateSamp(currentGenePool, popSize, samplingProb);
    #xCrossed = matrix(NA, nrow = popSize, ncol = geneLength)
    crossedSample <- matrix(NA, nrow = popSize, ncol = geneLength);
    #for(i in seq(1, popSize, by = 2))
    #  xCrossed[i:(i+1),] <- crossover(xSamp[i,], xSamp[i+1,], popSize, geneLength, crossRate)
    for(j in seq(1, popSize, by = 2)){
      #print(i)
      crossedSample[j:(j+1),] <- crossover(geneSample[j,], geneSample[j+1, ], geneLength, crossRate)
    }
    #
    #xMut = matrix(NA, nrow = popSize, ncol = geneLength)
    mutatedSample <- matrix(NA, nrow = popSize, ncol = geneLength)
    #for(i in seq(1, popSize, by = 2))
    #  xMut[i:(i+1),] <- mutation(xCrossed[i,], xCrossed[i+1,], mRate)  
    for (k in seq(1, popSize, by = 2)){
      #mutatedSample <- mutation(crossedSample[i,], crossedSample[i+1,], popSize, mRate)
      mutatedSample[k:(k+1),] <- mutation(crossedSample[k,], crossedSample[k+1,], mRate)
    }
    
    ### Here we would add the evaluation function ###
    # weights = AIC(  )
    currentGenePool <- mutatedSample
    samplingProb <- evalFunction(X, y, currentGenePool, popSize, type, family, criterion, criFun)[3,]
    avgCriterion <- rbind(avgCriterion, mean(evalFunction(X, y, currentGenePool, popSize, type, family, criterion, criFun)[1,]))
    #x = xMut # Update x-matrix with our new one!
    #print(x) # take out later
  }
  
  ##### After a fixed number of iterations, we return the best model #####
  #return(currentGenePool)
  final <- best(X, y, currentGenePool, popSize, type, criterion)
  #print(avgAIC)
  plot(avgCriterion, main='Average Criterion Values vs Iteration Number', xlab = "Iteration", ylab ="Average Criterion Value")
  ##### Print the best model #####
  return(final)
}

best <- function(X, y, pool, popSize, type, criterion, family = "gaussian", criFun = NULL){
  #print('In best')
  
  tmp <- evalFunction(X, y, pool, popSize, type, family, criterion, criFun)
  #print(result)
  final <- 0
  if(type == "lm"){
    #print('lm flow')
    index <- which(tmp[2,] == min(tmp[2,]), arr.ind = T)[1]
    #print(index)
    index2 <- which(pool[index,] != 0, arr.ind = T)
    #print(index2)
    final <- lm(y~X[,index2])
    #print('success')
  }
  else if (type == "glm"){
    #print('glm flow')
    index <- which(tmp[2,] == min(tmp[2,]), arr.ind = T)
    index2 <- which(pool[index,] != 0, arr.ind = T)[1]
    final <- glm(as.vector(y)~X[,index2],family)
    
  }
  ### TO be fixed here
  else{
    final <- 0
  }
  #print(summary(final))
  criFunBuilt <- eval(parse(text = criterion))
  criValue <- criFunBuilt(final)
  #print(paste("The resulting criterion is: ", criterion, criValue))
  return(final)
}




### test code
#result <- select(X, y, popSize = 19, max_iterations = 50, crossRate = 0.95, mRate = 0.0001)
v1 <- matrix(runif(200)*500,nrow = 200)
v2 <- matrix(runif(200)*10,nrow = 200)
error <- matrix(rnorm(200), nrow = 200)
n <- rep(200,20)
v3_22 <- sapply(n, runif)
v3_22 <- (v3_22)*500
X24 <- cbind(v1,v2,v3_22)
y1 <- 0.5*v1 + 30*v2 +error

