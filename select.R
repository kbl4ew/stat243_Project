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
X <- mtcars[,2:11]
y <- mtcars[,1]

select_test2 <- function(X = NULL, y = NULL, popSize = 200, criterion = "AIC", type = "lm", 
                   family = NA, criFun = NULL, max_iterations = 500, min_iterations = 50, 
                   crossRate = 1, mRate = NA, zeroToOneRatio = 2){
  
  ##### Defense coding #####
  X <- as.matrix(X)
  y <- as.vector(y)
  if(is.na(mRate)){
    mRate = 1/(dim(X)[1]);
  }
  
  if(is.null(evalFunction2)){
    stop("Please provide an evaluation function! Exiting from the function")
  }
  
  if(is.null(X)){
    stop("Please provide the predictors! Exiting from the function")
  }
  
  if(is.null(y)){
    stop("Please provide the dependent variable/outcome! Exiting from the function")
  }
  
  ##### Beginning of the generic algorithm #####
  geneLength <- dim(X)[2]
  
  ##### Initializing the first generation of individuals/ models
  initialPopulation <- popInitialize(popSize, geneLength, zeroToOneRatio)
  currentGenePool <- initialPopulation
  
  ### Calculating the sampling probabilities for the first generations of individuals/models
  #samplingProb <- evalFunction(currenGenePool, type, criterion, family, criFun)[3,]
  samplingProb <- evalFunction2(currentGenePool, popSize, type, family, criterion, criFun)[3,]
    
  ### While loop to handle convergence/ exceedance of min iteration/ capped by max iteration
  #iter = 0;
  ### Condition to be satisfied 
  ### if iter < min_iteration
  #while((iter <= min_iterations)&& !(iter >= max_iterations))
  
  for(i in 1:max_iterations){
    geneSample <- updateSamp(currentGenePool, popSize, samplingProb)
    crossedSample <- crossedParallel(geneSample, geneLength, crossRate, popSize)
    mutatedSample <- mutationParallel(crossedSample, mRate, popSize)
    
    currentGenePool <- mutatedSample
    samplingProb <- evalFunction2(currentGenePool, popSize, type, family, criterion, criFun)[3,]    
    
  }
  
  ##### After a fixed number of iterations, we return the best model #####
  return(currentGenePool)
  ##### Print the best model ##### 
}

#best <- function(currentGenePool, criterion, )

### test code
system.time({result_test <- select_test2(X, y, popSize = 50, max_iterations = 30, type = "lm", 
                                         crossRate = 0.95, mRate = 0.001)
             result})
