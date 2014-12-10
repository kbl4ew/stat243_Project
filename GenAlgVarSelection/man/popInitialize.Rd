\name{popInitialize}
\alias{popInitialize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
popInitialize(popSize = 0, geneLength = 0, zeroToOneRatio)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popSize}{
%%     ~~Describe \code{popSize} here~~
}
  \item{geneLength}{
%%     ~~Describe \code{geneLength} here~~
}
  \item{zeroToOneRatio}{
%%     ~~Describe \code{zeroToOneRatio} here~~
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
function (popSize = 0, geneLength = 0, zeroToOneRatio) 
{
    if (is.na(zeroToOneRatio)) {
        zeroToOneRation = 0
    }
    else {
    }
    pop <- matrix(nrow = popSize, ncol = geneLength)
    for (child in 1:popSize) {
        pop[child, ] = sample(c(rep(0, zeroToOneRatio), 1), geneLength, 
            replace = TRUE)
        while (sum(pop[child, ]) == 0) {
            pop[child, ] = sample(c(rep(0, zeroToOneRatio), 1), 
                geneLength, replace = TRUE)
        }
    }
    return(pop)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
