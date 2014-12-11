\name{select}
\alias{select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The variable selection algorithm
}

\description{
The Selection function uses the genetic algorithm to conduct variable selection for LMs and GLMs 
}
\usage{
select(X = NULL, y = NULL, popSize = 200, criterion = "AIC", type = "lm", 
          family = NA, criFun = NULL, max_iterations = 500, min_iterations = 50, 
          crossRate = NA, mRate = NA, zeroToOneRatio = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  User defined predictors, a matrix, of the linear model
}
  \item{y}{
  User defined response variable, a vector, of the linear model
}
  \item{popSize}{
  Number of models desired for each generation 
}
  \item{criterion}{
  The criterion used to evaluate individual models 
  (e.g. AIC, BIC)
}
  \item{type}{
  Type of regression of interest (i.e. lm or glm)
}
  \item{family}{
  A probability distribution from the exponential family to specialize the glm
}
  \item{criFun}{
  A function to evaluate the criterion when a built-in evaluation function is absent in R \cr
  The criFun is then passed to evalFunction to minimize via the Genetic Algorithm
}
  \item{max_iterations}{
  Maximum of iterations the Generic Algorithm carries out before stopping 
}
  \item{min_iterations}{
  Minimum of iterations the Generic Algorithm carries out before stopping
}
  \item{crossRate}{
  The probability that a sampled pair is crossovered
}
  \item{mRate}{
  The probability that each individual in the new generation is mutated
}
  \item{zeroToOneRatio}{
  The initial ratio of ommitted varaibles and included variables of the model
}
}
\details{
In the context of variable selection problems, each chromosome is encoded via binary variables. Given the nth chromosome, a 0 at the ith element indicates that the ith predictor is omitted while a 1 at the jth element indicates that the jth predictor is included in the nth chromosome/model.

The basic genetic algorithm is as follows:

Step 1: Initialize the first generation of models by generating a random population of n chromosomes;

Step 2: Evaluate the fitness/performance f(m) of each chromosome/model m in the first generation;

Step 3: Create a new population by repeating the following steps until the termination condition is met:

3a. Select n/2 pairs of parents chromosomes/models from the previous generation according to the fitness (fitter chromosomes have a greater chance to be selected)

3b. Crossover is carried out with a crossover probability. If crossover was performed, the 2 parent chromosomes are crossed over to produce 2 children. If no crossover was performed, 2 copies of the original chromosomes/models will be kept. 

3c. Mutation is carried out with a mutation probability. When mutation is carried out, the new offspring is mutated at each locus.

3d. Accept the children/offsprings as the next population

Step 4: If the terminal condition is satisfied, stop, and return the best solution in the current population. Otherwise, return to Step 2

Note: \cr
In order to avoid unnecessary computational cost, the authors have implemented a convergence check:
Once the average criterion value is no different that that of 10 iterations ago, the algorithm breaks.
}
\value{
Select() returns the best model (as evaluated by the user-defined criterion) after the algorithms reaches the terminal condition: i.e. if the uses inputs "lm" as the type, an "lm" object is returned; if the user inputs "glm" as the type, a "glm" object is returned.
}
\references{
Givens, G.H. and Hoeting, J.A. (2005). Computational Statistics. Wiley-Interscience \cr
}
\author{
Joy Hou, Kevin Li, Greta Olesen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Selecting the best predictors of MPG (using the built in mtcars dataset)
# Defining the predictors and response variables 
X <- mtcars[,2:11] 
y <- mtcars[,1] 
# Invoking select to return the best linear model 
result <- select(X, y, popSize = 19, max_iterations = 50, crossRate = 0.95, mRate = 0.0001)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
