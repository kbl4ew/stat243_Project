# The following lines are all for testing- take out later
## added variables: crossRate, mRate

set.seed(1)
x = matrix(sample(c(1,0), 30, replace = TRUE, prob=c(.7, .3)), nrow = 6)
popSize = dim(x)[1]
weights = (1:popSize)/sum(1:popSize) # would actually come from AIC/evaluation function
geneLength = dim(x)[2] 
crossRate = .9
mRate = .1

#############################################
#### 1st FUNCTION: SAMPLE FROM 'PARENTS' ####
#############################################

##### Input #####
# 1. x: 'Parent' matrix to current iteration
# 2. popSize: The population size in each generation
# 3. weights: Sampling weights produced by the AIC-evaluation function

##### Output #####
# 1. Updated matrix, ready for crossover and mutations

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

test1 = updateSamp(x, popSize, weights)
test1

############################################
######### 2nd FUNCTION: CROSSOVER ##########
############################################

##### Input #####
# 1. v1, v2: The function takes in one pair of models (for example: 1st row and 2nd row of matrix)
# 2. popSize: The population size in each generation
# 3. geneLength: The number of genes in the chromosome
# 4. crossRate: Rate at which crossover should be performed (defaults to 1)

##### Output #####
# 1. Post-crossover v1 and v2 in matrix of size 2 by popSize
# We do not want to any individual with all zeros hence we guarantee that at least one element is 1.

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
v1=(rep(0,5))
v2=(rep(1,5))
check = crossover(v1, v2, geneLength = length(v1)); check

############################################
######### 3rd FUNCTION: MUTATIONS ##########
############################################

##### Input #####
# 1. v1, v2: The function takes in one pair of models (for example: 1st row and 2nd row of matrix)
# 2. mRate: The rate of mutations
# 3. geneLength: the number of genes in the chromosome

##### Output #####
# 1. Post-mutated v1 and v2 in matrix of size 2 by popSize

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
    
    # old version
#    v1Loci = mLoci[mBool1] # mutation locations
#    v2Loci = mLoci[mBool2] # mutation locations 
#     v1Values = v1[v1Loci]   # value at the mutation locations
#     v1[v1Loci][v1Values == 1] <- 0
#     v1[v1Loci][v1Values == 0] <- 1
#     ## Change v2
#     v2Values = v2[v2Loci]   # value at the mutation locations
#     v2[v2Loci][v2Values == 1] <- 0
#     v2[v2Loci][v2Values == 0] <- 1
    
  }
}

v1 = rep(c(1,0), 3)
v2 = rep(1, 6)
check2 = mutation(v1, v2, mRate=1); check2 # mutates every time
check2 = mutation(v1, v2, mRate=0); check2 # never mutates


############################################
########## Implement Algorithm ############
############################################

# Make use of three functions above:
#   1) Implement the sampling-update function
#   2) Implement the crossover function
#   3) Implement the mutation function
#   4) Implement the evaluation/AIC function (not in this document)
#   5) Go back to Step 1 with updated matrix and new weights

for(i in 1:3){ 
  # really we will have predetermined # of iterations
  xSamp <- updateSamp(x, popSize, weights)
  
  xCrossed = matrix(NA, nrow = popSize, ncol = geneLength)
  for(i in seq(1, popSize, by = 2))
    xCrossed[i:(i+1),] <- crossover(xSamp[i,], xSamp[i+1,], popSize, geneLength)
  
  xMut = matrix(NA, nrow = popSize, ncol = geneLength)
  for(i in seq(1, popSize, by = 2))
    xMut[i:(i+1),] <- mutation(xCrossed[i,], xCrossed[i+1,], mRate)  
  
  ### Here we would add the evaluation function ###
  # weights = AIC(  )
  
  x = xMut # Update x-matrix with our new one!
  print(x) # take out later
}


### Notes: The way I envision making the code faster is to use parallel processing 
#   in the 2 individual for-loops that implement 'crossover' and 'mutations'
