##### Population OOP Code #####
##### Description: #####
### Here we construct the basic S4 Class: Population ###

setClass("Population", 
         representation(encoding = "character", popSize = "numeric", geneLength = "numeric", 
                        genePool = "matrix", fit = "vector", generation = "numeric"), 
         prototype(encoding = "binary", popSize = 0, geneLength = 0, genePool = matrix(nrow = 0, ncol = 0),
                   fit = c(), generation = 0)
)

##### Constructor #####


population <- function(encoding = "binary", Size = 200, Length = 200, initializeFunction = NA,...){
  ### defense mechanism ###
  ### To be done ###
  
  ### Assignment ###
  newObject <- new("Population")
  
  newObject@encoding = "binary"
  newObject@popSize = Size
  newObject@geneLength = Length
  newObject@genePool = initializeFunction(Size, Length)
  newObject@fit = rep(0, Size)  
  newObject@generation = 0
}




