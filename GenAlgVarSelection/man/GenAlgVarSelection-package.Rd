\name{GenAlgVarSelection-package}
\alias{GenAlgVarSelection}
\docType{package}
\title{
GenAlgVarSelecton
}
\description{
This package utilizes the genetic algorithm to conduct variable selection for Linear Models and Generalized Linear Models
}
\details{
\tabular{ll}{
Package: \tab GenAlgVarSelection\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-09\cr
Depends: \tab parallel, doParallel, foreach, iterators \cr
}
}
The GenAlgVarSelection Package provides a variable selection algorithm for Linear Models and Generalized Linear Models. The main function of the package is select().

\author{
Joy Hou, Kevin Li, Greta Olesen \cr
Maintainer: \cr
Joy Hou <juehou@berkeley.edu>, Kevin Li <kbl4ew@berkeley.edu>, Greta Olesen <greta.olesen@berkeley.edu>\cr
}

References:
\references{
Givens, G.H. and Hoeting, J.A. (2005). Computational Statistics. Wiley-Interscience \cr
}

Examples: \cr
\examples{
# Selecting the best predictors of MPG (using the built in mtcars dataset) 
# Defining the predictors and response variables 
X <- mtcars[,2:11]
y <- mtcars[,1] 
# Invoking select to return the best linear model 
result <- select(X, y, popSize = 19, max_iterations = 50, crossRate = 0.95, mRate = 0.0001)
}

