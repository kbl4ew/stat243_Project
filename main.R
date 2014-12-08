############################ Project: Main Code #################################

# Make sure to run these files first:
#   1) exploratory_code.R
#   2) parallelization.R
#   3) evaluationFunction.R
#   4) popInitualize.R

# parameters of the main function:
# ...
# ...
# model: either 1 (for lm) or 2(for glm)
# eval: evaluation criterion

########## ALGORITHM ###############
genAlg <- function(covariates, outcome, x, popSize, geneLength, crossRate, mRate, model,
                   family = gaussian, criterion = "AIC", criFun){
  if ((model ! = 1)&&(model != 2) ){
    stop("Model needs to be 1(lm) or 2(glm).")
  }
  
  
  for(i in 1:3){      # really we will have predetermined # of iterations
    
    if(i == 1){ # ADD THIS!!!!!!!!!!
      #x = popInitialize(    )
      #weights = evalFunction(genePool = xMut, covariates, outcome)[3,]
    }
    
    xSamp <- updateSamp(x, popSize, weights)
    
    xCrossed = crossedParallel(xSamp, geneLength, crossRate)
    xMut = mutationParallel(xCrossed, mRate)
    
    x = xMut # Update x-matrix with our new one!
    
    if(model = 1){
      weights <- evalLm(genePool = xMut, covariates, outcome, criterion, criFun)[3,]
    } else {
      weights <- evalGlm(genePool = xMut, covariates, outcome, family, criterion, criFun)[3,]
    }
    
    
    print(x) # take out later
  }
}
genAlg(X, y, x, popSize, geneLength, crossRate, mRate)

