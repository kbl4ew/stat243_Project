library(parallel) 
library(doParallel)
library(foreach)
library(iterators)
nCores <- 4  
registerDoParallel(nCores) 

crossedParallel = function(x, geneLength, crossRate){
  foreach(i=seq(1, popSize, by = 2), .combine = rbind) %dopar% {
    crossedPair <- crossover(x[i,], x[i+1,], geneLength, crossRate)
    return(crossedPair)
  } 
}

## tests
v1=(rep(0,5)); v2=(rep(1,5))
x2 = rbind(v1,v2,v1,v2,v1,v2)
crossedParallel(x, geneLength, crossRate)


mutationParallel = function(x, mRate){
  foreach(i=seq(1, popSize, by = 2), .combine = rbind) %dopar% {
    mutatedPair <- mutation(x[i,], x[i+1,], mRate)
    return(mutatedPair)
  } 
}

## TESTS
v1 = rep(c(1,0), 3); v2 = rep(1, 6)
x2 = rbind(v1,v2,v1,v2,v1,v2)
check2 = mutationParallel(x2, mRate=1); check2 # mutates every time
check2 = mutationParallel(x2, mRate=0); check2 # never mutates

########## ALGORITHM ###############

genAlg <- function(x, popSize, geneLength, crossRate, mRate){
  for(i in 1:3){      # really we will have predetermined # of iterations
    ### Here we would add the evaluation function ###
    # weights = AIC(x  )
    
    xSamp <- updateSamp(x, popSize, weights)
    
    xCrossed = crossedParallel(xSamp, geneLength, crossRate)
    xMut = mutationParallel(xCrossed, mRate)
    
    #add weights here?
    x = xMut # Update x-matrix with our new one!
    print(x) # take out later
  }
}
genAlg(x, popSize, geneLength, crossRate, mRate)
