\name{crossover}
\alias{crossover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
crossover(v1, v2, geneLength, crossRate = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v1}{
%%     ~~Describe \code{v1} here~~
}
  \item{v2}{
%%     ~~Describe \code{v2} here~~
}
  \item{geneLength}{
%%     ~~Describe \code{geneLength} here~~
}
  \item{crossRate}{
%%     ~~Describe \code{crossRate} here~~
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
function (v1, v2, geneLength, crossRate = 1) 
{
    crossBool = sample(c(TRUE, FALSE), 1, prob = c(crossRate, 
        1 - crossRate))
    if (crossBool) {
        cut = sample(geneLength - 1, 1)
        new1 = c(v1[1:cut], v2[(cut + 1):geneLength])
        new2 = c(v2[1:cut], v1[(cut + 1):geneLength])
        if (sum(new1) == 0 | sum(new2) == 0) 
            return(rbind(v1, v2))
        else return(rbind(new1, new2))
    }
    else return(rbind(v1, v2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
