\name{singleEval}
\alias{singleEval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
singleEval(singleGene, X, y, type, criterion, criFun, family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{singleGene}{
%%     ~~Describe \code{singleGene} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
}
  \item{criFun}{
%%     ~~Describe \code{criFun} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
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
function (singleGene, X, y, type, criterion, criFun, family) 
{
    if (type == "lm") 
        fit <- lm(y ~ as.matrix(X[, which(singleGene != 0)]))
    if (type == "glm") 
        fit <- glm(as.vector(y) ~ as.matrix(X[, which(singleGene != 
            0)]), family)
    if (is.null(criFun)) {
        criFunBuilt <- eval(parse(text = criterion))
        criValue <- criFunBuilt(fit)
    }
    else criValue = criFun(fit)
    return(criValue)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
