\name{select}
\alias{select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
The Selection function uses the genetic algorithm to conduct variable selection for LMs and GLMs 
}
\usage{
select(X = NULL, y = NULL, popSize = 200, criterion = "AIC", type = "lm", family = NA, criFun = NULL, max_iterations = 500, min_iterations = 50, crossRate = NA, mRate = NA, zeroToOneRatio = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  User defined predictors of the linear model
}
  \item{y}{
  User defined response variable of the linear model
}
  \item{popSize}{
  Number of models desired for each generation 
}
  \item{criterion}{
  The criterion used to evaluate individual models
}
  \item{type}{
  Type of regression of interest (i.e. lm or glm)
}
  \item{family}{
  A probability distribution from the exponential family to specialize the glm
}
  \item{criFun}{
  A function to evaluate the criterion when 
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
  The initial ratio of ommitted varaibles and included variables
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X = NULL, y = NULL, popSize = 200, criterion = "AIC", 
    type = "lm", family = NA, criFun = NULL, max_iterations = 500, 
    min_iterations = 50, crossRate = NA, mRate = NA, zeroToOneRatio = 2) 
{
    X <- as.matrix(X)
    y <- as.vector(y)
    if (is.na(mRate)) {
        mRate = 1/(dim(X)[1])
    }
    if ((popSize\%\%2) != 0) {
        print("Warning: The number of models has been incremented to the nearest even number")
        warning("The number of models has been incremented to the nearest even number")
        popSize <- popSize + 1
    }
    if (is.null(X)) {
        stop("Please provide the predictors! Exiting from the function")
    }
    if (is.null(y)) {
        stop("Please provide the independent variable/outcome! Exiting from the function")
    }
    geneLength <- dim(X)[2]
    initialPopulation <- popInitialize(popSize, geneLength, zeroToOneRatio)
    currentGenePool <- initialPopulation
    samplingProb <- evalFunction(currentGenePool, popSize, type, 
        family, criterion, criFun)[3, ]
    avgCriterion <- mean(evalFunction(currentGenePool, popSize, 
        type, family, criterion, criFun)[1, ])
    for (i in 1:max_iterations) {
        geneSample <- updateSamp(currentGenePool, popSize, samplingProb)
        crossedSample <- matrix(NA, nrow = popSize, ncol = geneLength)
        for (i in seq(1, popSize, by = 2)) {
            crossedSample[i:(i + 1), ] <- crossover(geneSample[i, 
                ], geneSample[i + 1, ], geneLength, crossRate)
        }
        mutatedSample <- matrix(NA, nrow = popSize, ncol = geneLength)
        for (i in seq(1, popSize, by = 2)) {
            mutatedSample[i:(i + 1), ] <- mutation(crossedSample[i, 
                ], crossedSample[i + 1, ], mRate)
        }
        currentGenePool <- mutatedSample
        samplingProb <- evalFunction(currentGenePool, popSize, 
            type, family, criterion, criFun)[3, ]
        avgCriterion <- rbind(avgCriterion, mean(evalFunction(currentGenePool, 
            popSize, type, family, criterion, criFun)[1, ]))
    }
    final <- best(currentGenePool, popSize, type, criterion)
    plot(avgCriterion)
    return(final)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
