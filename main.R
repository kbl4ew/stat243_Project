############################ Project: Main Code #################################

# Make sure to run these files first:
#   1) exploratory_code.R
#   2) parallelization.R
#   3) evaluationFunction.R
#   4) popInitualize.R

########## ALGORITHM ###############
genAlg <- function(covariates, outcome, x, popSize, geneLength, cross_rate, m_rate){
  for(i in 1:3){      # really we will have predetermined # of iterations
    
    if(i == 1){ # ADD THIS!!!!!!!!!!
      #x = popInitialize(    )
      #weights = evalFunction(genePool = x_mut, covariates, outcome)[3,]
    }
    
    x_samp <- update_samp(x, popSize, weights)
    
    x_crossed = crossed_parallel(x_samp, geneLength, cross_rate)
    x_mut = mutation_parallel(x_crossed, m_rate)
    
    x = x_mut # Update x-matrix with our new one!
    weights <- evalFunction(genePool = x_mut, covariates, outcome)[3,]
    
    print(x) # take out later
  }
}
genAlg(X, y, x, popSize, geneLength, cross_rate, m_rate)

