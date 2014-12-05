############## Code to randomize the first generation for generic algorithm #####################

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
popInitialize <- function(popSize = 0, geneLength = 0, zeroToOneRatio){
  
  if(is.na(zeroToOneRatio)){
    zeroToOneRation = 0;    
  }
  else{
    print('zeroToOneRatio accepted')
  }
  pop <- matrix(nrow = popSize, ncol = geneLength);

  ##### Randomly initialize the first generation #####
  for (child in 1:popSize){
    pop[child, ] = sample(c(rep(0, zeroToOneRatio), 1), geneLength, replace = TRUE);
    #print(child)
    #print(pop[child, ])
    while(sum(pop[child,]) == 0){
      pop[child, ] = sample(c(rep(0, zeroToOneRatio), 1), geneLength, replace = TRUE);   
      #print('in the loop')
      #print(pop[child, ])
    }
  }
  return(pop)
}

##### Test Case #####
tmp <- popInitialize(6, geneLength = 4, zeroToOneRatio = 2 )
