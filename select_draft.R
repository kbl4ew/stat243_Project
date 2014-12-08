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

select <- function(X = NULL, y = NULL, popSize = 200, criterion = "AIC", type = "lm", family = NA, evalFunction = NULL,max_iterations = 500, min_iterations = 50, mutation_prob = NA,  zeroToOneRatio = 10){
  ##### Defense coding #####
  X <- matrix(X)
  y <- as.vector(y)
  if(is.na(mutation_Prob)){
    mutation_Prob = 1/(dim(X)[1]);
  }
  
  if(is.null(evalFunction)){
    stop("Please provide an evaluation function! Exiting from the function")
  }
  
  if(is.null(X)){
    stop("Please provide the predictors! Exiting from the function")
  }
  
  if(is.null(y)){
    stop("Please provide the independent variable/outcome! Exiting from the function")
  }
  
  ##### Beginning of the generic algorithm #####
  geneLength <- dim(X)[2]
  ##### Initializing the first generation of individuals/ models
  initial_population <- popInitialize(popSize = 0, geneLength, zeroToOneRatio)
  
  ### Calculating the sampling probabilities for the first generations of individuals/models
  samplingProb <- evalFunction(type, criterion, family, evalFunction)
}




















