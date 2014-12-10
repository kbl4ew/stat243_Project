\name{mutation}
\alias{mutation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mutation(v1, v2, mRate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v1}{
%%     ~~Describe \code{v1} here~~
}
  \item{v2}{
%%     ~~Describe \code{v2} here~~
}
  \item{mRate}{
%%     ~~Describe \code{mRate} here~~
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
function (v1, v2, mRate) 
{
    mLoci = which(v1 == v2)
    len = length(mLoci)
    mBool1 = sample(c(TRUE, FALSE), len, replace = TRUE, prob = c(mRate, 
        1 - mRate))
    mBool2 = sample(c(TRUE, FALSE), len, replace = TRUE, prob = c(mRate, 
        1 - mRate))
    if (sum(mBool1) == 0 & sum(mBool2) == 0) 
        return(rbind(v1, v2))
    else {
        v1Copy = v1
        v2Copy = v2
        v1Copy[mLoci][mBool1] <- as.numeric(!v1Copy[mLoci][mBool1])
        v2Copy[mLoci][mBool2] <- as.numeric(!v2Copy[mLoci][mBool2])
        if (sum(v1Copy) == 0) 
            v1Copy <- v1
        if (sum(v2Copy) == 0) 
            v2Copy <- v2
        return(rbind(v1Copy, v2Copy))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
