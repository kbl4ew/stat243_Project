library(parallel) 
library(doParallel)
library(foreach)
library(iterators)
nCores <- 4  
registerDoParallel(nCores) 

crossed_parallel = function(x, geneLength, cross_rate){
  foreach(i=seq(1, popSize, by = 2), .combine = rbind) %dopar% {
    crossed_pair <- crossover(x[i,], x[i+1,], geneLength, cross_rate)
    return(crossed_pair)
  } 
}

## tests
v1=(rep(0,5)); v2=(rep(1,5))
x2 = rbind(v1,v2,v1,v2,v1,v2)
crossed_parallel(x, geneLength, cross_rate)


mutation_parallel = function(x, m_rate){
  foreach(i=seq(1, popSize, by = 2), .combine = rbind) %dopar% {
    mutated_pair <- mutation(x[i,], x[i+1,], m_rate)
    return(mutated_pair)
  } 
}

## TESTS
v1 = rep(c(1,0), 3); v2 = rep(1, 6)
x2 = rbind(v1,v2,v1,v2,v1,v2)
check2 = mutation_parallel(x2, m_rate=1); check2 # mutates every time
check2 = mutation_parallel(x2, m_rate=0); check2 # never mutates

########## ALGORITHM ###############

genAlg <- function(x, popSize, geneLength, cross_rate, m_rate){
  for(i in 1:3){      # really we will have predetermined # of iterations
    ### Here we would add the evaluation function ###
    # weights = AIC(x  )
    
    x_samp <- update_samp(x, popSize, weights)
    
    x_crossed = crossed_parallel(x_samp, geneLength, cross_rate)
    x_mut = mutation_parallel(x_crossed, m_rate)
    
    #add weights here?
    x = x_mut # Update x-matrix with our new one!
    print(x) # take out later
  }
}
genAlg(x, popSize, geneLength, cross_rate, m_rate)
