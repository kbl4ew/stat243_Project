\name{GenAlgVarSelection-package}
\alias{GenAlgVarSelection-package}
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
The GenAlgVarSelection Package provides a variable selection algorithm for Linear Models and Generalized Linear Models. The main function of the package is select().

\author{
Joy Hou, Kevin Li, Greta Olesen
Maintainer: Kevin Li <kbl4ew@berkeley.edu> \cr
}
References:
\references{
Givens, G.H. and Hoeting, J.A. (2005). Computational Statistics. Wiley-Interscience \cr
}

Examples: \cr
\examples{
# Selecting the best predictors of MPG (using the built in mtcars dataset) \cr
# Defining the predictors and response variables \cr
X <- mtcars[,2:11] \cr
y <- mtcars[,1] \cr
# Invoking select to return the best linear model \cr
result <- select(X, y, popSize = 19, max_iterations = 50, crossRate = 0.95, mRate = 0.0001)
}

