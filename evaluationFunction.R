##### Evaluation Function #####
##### Descritption ######
### Here we create an evalFunction that allows us evaluate the fitness of the current generation. ###

##### Input #####
### genePool: a matrix holding the current individuals of size population size by 
### X: the independent variables (training data)
### y: the dependent varaibles (training data)

##### Output #####
### result: a matrix which holds the AIC & rank (size: 2 by populationSize)
### First row of result: AIC
### Second row of result rank
#library(data.table)

evalFunction <- function(genePool, X, y){
  popSize <- dim(genePool)[1] 
  geneLength <- dim(genePool)[2]
  result <- matrix(nrow = 1, ncol = popSize)
  rownames(result) <- c("AIC")
  ### a for loop to iterate across each individual to calculate 
  for (i in 1:popSize){
    fit <- lm(y~as.matrix(X[,which(genePool[i,] != 0, arr.ind = T)]))
    result[1, i] <- AIC(fit)
    #result[1, i] <- -i
  }
#   DT <- as.data.table(t(result))
#   DT[, AIC_Rank := rank(AIC, ties.method = "first")]
#   return(t(DT))
  
  obj <- rbind(result, rank(result), rank(result)/sum(1:popSize))
  row.names(obj) <- c("AIC", "ranks", "sampling_probs")
  return(obj)
}


evalFunction()
dim(mtcars)

evalFunction(test, X, y)
#rank(AIC_vals) / sum(1:length(AIC_vals)
