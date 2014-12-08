data <- mtcars[,2:10]
y <- mtcars[,1]
popSize = 20 # this can be anything, as long as it's even!
geneLength = dim(data)[2] 
crossRate = .9
mRate = .1
weights <- rep(1,popSize)/popSize

itrNum = 100

###test popInitialize()
genePool <- popInitialize(popSize, geneLength, zeroToOneRatio = 2)

for(i in 1:itrNum){ 
  # really we will have predetermined # of iterations
  xSamp <- updateSamp(genePool, popSize, weights)
  
  xCrossed = matrix(NA, nrow = popSize, ncol = geneLength)
  for(i in seq(1, popSize, by = 2))
    xCrossed[i:(i+1),] <- crossover(xSamp[i,], xSamp[i+1,],geneLength,crossRate)
  
  xMut = matrix(NA, nrow = popSize, ncol = geneLength)
  for(i in seq(1, popSize, by = 2))
    xMut[i:(i+1),] <- mutation(xCrossed[i,], xCrossed[i+1,], mRate)  
  
  ### Here we would add the evaluation function ###
  if(model = 1){
    weights <- evalLm(genePool = xMut, covariates, outcome, criterion, criFun)[3,]
  } else {
    weights <- evalGlm(genePool = xMut, covariates, outcome, family, criterion, criFun)[3,]
  }
  
  x = xMut # Update x-matrix with our new one!
  print(x) # take out later
}
