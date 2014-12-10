\name{best}
\alias{best}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
best(pool, popSize, type, criterion, family = NA, criFun = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pool}{
%%     ~~Describe \code{pool} here~~
}
  \item{popSize}{
%%     ~~Describe \code{popSize} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
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
function (pool, popSize, type, criterion, family = NA, criFun = NULL) 
{
    tmp <- evalFunction(pool, popSize, type, family, criterion, 
        criFun)
    final <- 0
    if (type == "lm") {
        index <- which(tmp[2, ] == min(tmp[2, ]), arr.ind = T)[1]
        index2 <- which(pool[index, ] != 0, arr.ind = T)
        final <- lm(y ~ as.matrix(X[, index2]))
    }
    else if (type == "glm") {
        index <- which(tmp[2, ] == min(tmp[2, ]), arr.ind = T)
        index2 <- which(pool[index, ] != 0, arr.ind = T)[1]
        final <- glm(as.vector(y) ~ as.matrix(X[, index2]), family)
    }
    else {
        final <- 0
    }
    print(final)
    print(paste("The resulting criterion is: ", criterion, AIC(final)))
    return(final)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
