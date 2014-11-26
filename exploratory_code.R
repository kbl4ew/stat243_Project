# The following lines are all for testing- take out later
## added variables: cross_rate, m_rate

set.seed(1)
x = matrix(sample(c(1,0), 24, replace = TRUE, prob=c(.7, .3)), nrow = 6)
popSize = dim(x)[1]
weights = (1:popSize)/sum(1:popSize) # would actually come from AIC/evaluation function
geneLength = dim(x)[2] 
cross_rate = .9
m_rate = .1
v1 = x[1,]
v2 = x[2,]

#############################################
#### 1st FUNCTION: SAMPLE FROM 'PARENTS' ####
#############################################

##### Input #####
# 1. x: 'Parent' matrix to current iteration
# 2. popSize: The population size in each generation
# 3. weights: Sampling weights produced by the AIC-evaluation function

##### Output #####
# 1. Updated matrix, ready for crossover and mutations

update_samp <- function(x, popSize, weights){
  pairs = matrix(NA, ncol = popSize/2, nrow = 2) # form sampled pairs
  
  set.seed(0) # take out later
  for(i in 1:popSize/2){
    pairs[,i] = sample(1:popSize, 2, prob = weights)
  }
  pairs # so 2 and 6 are paired up; 4 and 5 paired up; 6,5 paired
  x_samp = x[as.vector(pairs), ]
  return(x_samp)
}

test1 = update_samp(x, popSize, weights)
test1

############################################
######### 2nd FUNCTION: CROSSOVER ##########
############################################

##### Input #####
# 1. v1, v2: The function takes in one pair of models (for example: 1st row and 2nd row of matrix)
# 2. popSize: The population size in each generation
# 3. geneLength: The number of genes in the chromosome
# 4. cross_rate: Rate at which crossover should be performed (defaults to 1)

##### Output #####
# 1. Post-crossover v1 and v2 in matrix of size 2 by popSize
# We do not want to any individual with all zeros hence we guarantee that at least one element is 1.

crossover <- function(v1, v2, popSize, geneLength, cross_rate=1){
  cross_bool = sample(c(TRUE, FALSE), 1, prob = c(cross_rate, 1-cross_rate))
  
  if(cross_bool){
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

test=matrix(NA, nrow=2, ncol = 4)
test[1:2,] <- crossover(v1, v2, popSize, geneLength, cross_rate); test

############################################
######### 3rd FUNCTION: MUTATIONS ##########
############################################

##### Input #####
# 1. v1, v2: The function takes in one pair of models (for example: 1st row and 2nd row of matrix)
# 2. m_rate: The rate of mutations
# 3. geneLength: the number of genes in the chromosome

##### Output #####
# 1. Post-mutated v1 and v2 in matrix of size 2 by popSize
# TO DO: guarantee that at least one element is 1.

mutation <- function(v1, v2, m_rate){
  m_loci = which(v1==v2)
  len = length(m_loci)
  
  # T/F: mutate or not
  m_bool1 = sample(c(TRUE, FALSE), len, replace = TRUE, prob = c(m_rate, 1-m_rate))
  m_bool2 = sample(c(TRUE, FALSE), len, replace = TRUE, prob = c(m_rate, 1-m_rate))
            
  if(sum(m_bool1) == 0 & sum(m_bool2) == 0) 
    return(rbind(v1,v2)) # return v1 v2 and dont mutate
  
  else{
    ## Change v1:
    v1_loci = m_loci[m_bool1] # mutation locations
    v1_values = v1[v1_loci]   # value at the mutation locations
    v1[v1_loci][v1_values == 1] <- 0
    v1[v1_loci][v1_values == 0] <- 1
    
    ## Change v2
    v2_loci = m_loci[m_bool2] # mutation locations
    v2_values = v2[v2_loci]   # value at the mutation locations
    v2[v2_loci][v2_values == 1] <- 0
    v2[v2_loci][v2_values == 0] <- 1
    
    return(rbind(v1,v2))
  }
}

test2=matrix(NA, nrow=2, ncol = 4)
test2[1:2,] <- mutation(v1, v2, m_rate); test2


############################################
########## Implement Algorithm ############
############################################

# Make use of three functions above:
#   1) Implement the sampling-update function
#   2) Implement the crossover function
#   3) Implement the mutation function
#   4) Implement the evaluation/AIC function (not in this document)
#   5) Go back to Step 1 with updated matrix and new weights

for(i in 1:3){ # really we will have predetermined # of iterations
  x_samp <- update_samp(x, popSize, weights)
  
  x_crossed = matrix(NA, nrow = popSize, ncol = geneLength)
  for(i in seq(1, popSize, by = 2))
    x_crossed[i:(i+1),] <- crossover(x_samp[i,], x_samp[i+1,], popSize, geneLength, cross_rate)
  
  x_mut = matrix(NA, nrow = popSize, ncol = geneLength)
  for(i in seq(1, popSize, by = 2))
    x_mut[i:(i+1),] <- mutation(x_crossed[i,], x_crossed[i+1,], m_rate)  
  
  ### Here we would add the evaluation function ###
  # weights = AIC(  )
  
  x = x_mut # Update x-matrix with our new one!
  print(x) # take out later
}


### Notes: The way I envision making the code faster is to use parallel processing 
#   in the 2 individual for-loops that implement 'crossover' and 'mutations'
