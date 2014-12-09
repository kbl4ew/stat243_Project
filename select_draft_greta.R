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
X <- mtcars[,2:11]
y <- mtcars[,1]
select_test <- function(X = NULL, y = NULL, popSize = 200, criterion = "AIC", type = "lm", family = NA, criFun = NULL, max_iterations = 500, min_iterations = 50, crossRate = NA, mRate = NA, zeroToOneRatio = 2){
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
    stop("Please provide the independent variable/outcome! Exiting from the function")
  }
  
  ##### Beginning of the generic algorithm #####
  geneLength <- dim(X)[2]
  ##### Initializing the first generation of individuals/ models
  initialPopulation <- popInitialize(popSize, geneLength, zeroToOneRatio)
  currentGenePool <- initialPopulation
  ### Calculating the sampling probabilities for the first generations of individuals/models
  #samplingProb <- evalFunction2(currenGenePool, type, criterion, family, criFun)[3,]
    
  samplingProb <- evalFunction2(currentGenePool, popSize, type, family, criterion, criFun)[3,]
  ### While loop to handle convergence/ exceedance of min iteration/ capped by max iteration
  #iter = 0;
  ### Condition to be satisfied 
  ### if iter < min_iteration
  #while((iter <= min_iterations)&& !(iter >= max_iterations))
  for(i in 1:max_iterations){
    # really we will have predetermined # of iterations
    #xSamp <- updateSamp(x, popSize, weights)
    geneSample <- updateSamp(currentGenePool, popSize, samplingProb)
    #xCrossed = matrix(NA, nrow = popSize, ncol = geneLength)
    crossedSample <- matrix(NA, nrow = popSize, ncol = geneLength)
    #for(i in seq(1, popSize, by = 2))
    #  xCrossed[i:(i+1),] <- crossover(xSamp[i,], xSamp[i+1,], popSize, geneLength, crossRate)
    for(i in seq(1, popSize, by = 2)){
      crossedSample[i:(i+1),] <- crossover(geneSample[i,], geneSample[i+1, ], geneLength, crossRate)
    }
    #
    #xMut = matrix(NA, nrow = popSize, ncol = geneLength)
    mutatedSample <- matrix(NA, nrow = popSize, ncol = geneLength)
    #for(i in seq(1, popSize, by = 2))
    #  xMut[i:(i+1),] <- mutation(xCrossed[i,], xCrossed[i+1,], mRate)  
    for (i in seq(1, popSize, by = 2)){
      #mutatedSample <- mutation(crossedSample[i,], crossedSample[i+1,], popSize, mRate)
      mutatedSample[i:(i+1),] <- mutation(crossedSample[i,], crossedSample[i+1,], mRate)
    }

    ### Here we would add the evaluation function ###
    # weights = AIC(  )
    
    currentGenePool <- mutatedSample
    
    ### HERE IS WHERE IT IS DIFFERENT
    samplingProb <- evalFunction2(currentGenePool, popSize, type, family, criterion, criFun)[3,]
    #x = xMut # Update x-matrix with our new one!
  }
  
  ##### After a fixed number of iterations, we return the best model #####
  return(currentGenePool)
  
  ##### Print the best model #####
  
}
system.time({result <- select_test(X, y, popSize = 200, max_iterations = 100, crossRate = 0.95, mRate = 0.001)
             result})

# select is in the 
system.time({result <- select(X, y, popSize = 200, max_iterations = 100, crossRate = 0.95, mRate = 0.001)
             result})
result <- select(X, y, popSize = 20, max_iterations = 50, crossRate = 0.95, mRate = 0.001)


#best <- function(currentGenePool, criterion, )




### test code
result <- select(X, y, popSize = 20, max_iterations = 50, crossRate = 0.95, mRate = 0.001)






X
y




