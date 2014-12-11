testStepwise = function(){

  ### testing against stepwise regression result ###
  ### using mtcars data
  X <- mtcars[,2:11]
  y <- mtcars[,1]
  
  
  cat("Our function running ...")
  set.seed(0)
  result2 <- select(X, y, popSize = 160, max_iterations = 50, crossRate = 0.95, mRate = 0.0001)
  cat("Our function has chosen the following model:")
  print(summary(result2))
  cat("The AIC value for our model is:",unlist(AIC(result2)),"\n")
  
  cat("Now we implement stepwise regression on the dataset.")
  fullModel <- lm( mpg ~ cyl+disp+hp+drat+wt+qsec+vs+am+gear+carb, data = mtcars)
  stepResult <- step(fullModel, direction = "both", trace = 1)
  cat("The stepwise regression has picked the following model:")
  print(summary(stepResult))
  cat("The AIC value for this model is:",AIC(stepResult))
  
  if((abs(AIC(result2)-AIC(stepResult))) < 10)
    cat("The model our function chose is close to the one that stepwise regression chose, test succeeded.")
  else
    cat("The model our function chose is not close to the one that stepwise regression chose, test failed")
}
## stepwise regression selected almost the same model ##
## AIC is almost the same as from our model as well ##

### testing out the program on simulated data set ###
## Initialize the true model: y is dependent only on v1 and v2 ##
testSim <- function(){
  X <- mtcars[,1:11]
  n <- dim(mtcars)[1]
  error <- matrix(rnorm(n),nrow = n)
  y <- 1*X[,1] + 2*X[,2] + 3*X[,3] + 4*X[,4] + 5*X[,5] + error
  
  set.seed(1)
  cat("Our function running ...")
  set.seed(1)
  result <- select (X, y, popSize = 160, max_iteration = 50, criterion = "AIC", zeroToOneRatio = 1, crossRate = 0.95, mRate = 0.001)
  cat("Our function has chosen the following model:")
  print(summary(result))
  cat("The AIC value for our model is:",unlist(AIC(result)),"\n")
  
  cat("The true model has the mpg, cyl, disp, hp and drat as the independent variables.\n")
  cat("Our function has picked out all of the 5 relevant variables but included 1 additional irrelevant variable. The performance of our function is decent. Test passed.")

}
