\name{evalFunction}
\alias{evalFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evalFunction(currentGenePool, popSize, type = "lm", family = "gaussian", criterion = "AIC", criFun = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{currentGenePool}{
%%     ~~Describe \code{currentGenePool} here~~
}
  \item{popSize}{
%%     ~~Describe \code{popSize} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
}
  \item{criFun}{
%%     ~~Describe \code{criFun} here~~
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
%%  ~~who you are~~
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
function (currentGenePool, popSize, type = "lm", family = "gaussian", 
    criterion = "AIC", criFun = NULL) 
{
    if ((popSize\%\%2) != 0) {
        message("Warning: The size of the population has been rounded to the largest even numer")
        popSize <- popSize + 1
    }
    if (criterion != "AIC" & criterion != "BIC") 
        stop(paste(criterion, "is not a valid criterion. Please use AIC or BIC."))
    if (!is.null(criFun) & !is.function(criFun)) 
        stop("criFun input is not a function.")
    if (type != "lm" & type != "glm") 
        stop("Regression must be of type 'lm' or 'glm'")
    geneLength <- dim(currentGenePool)[2]
    result <- rep(NA, popSize)
    result <- foreach(i = 1:popSize, .combine = c) \%dopar\% {
        criValue <- singleEval(currentGenePool[i, ], X, y, type, 
            criterion, criFun, family)
        return(criValue)
    }
    obj <- rbind(result, rank(result), rank(-result)/sum(1:popSize))
    row.names(obj) <- c(criterion, "ranks", "samplingProbs")
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
