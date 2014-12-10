\name{result}
\alias{result}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(result)}
\format{
  The format is:
List of 12
 $ coefficients : Named num [1:4] 12.943 -4.618 1.049 0.891
  ..- attr(*, "names")= chr [1:4] "(Intercept)" "as.matrix(X[, index2])wt" "as.matrix(X[, index2])qsec" "as.matrix(X[, index2])gear"
 $ residuals    : Named num [1:32] -0.6709 -0.0806 -2.5108 0.2432 1.1199 ...
  ..- attr(*, "names")= chr [1:32] "1" "2" "3" "4" ...
 $ effects      : Named num [1:32] -113.65 -29.116 -9.103 2.735 0.952 ...
  ..- attr(*, "names")= chr [1:32] "(Intercept)" "as.matrix(X[, index2])wt" "as.matrix(X[, index2])qsec" "as.matrix(X[, index2])gear" ...
 $ rank         : int 4
 $ fitted.values: Named num [1:32] 21.7 21.1 25.3 21.2 17.6 ...
  ..- attr(*, "names")= chr [1:32] "1" "2" "3" "4" ...
 $ assign       : int [1:4] 0 1 1 1
 $ qr           :List of 5
  ..$ qr   : num [1:32, 1:4] -5.657 0.177 0.177 0.177 0.177 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:32] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:4] "(Intercept)" "as.matrix(X[, index2])wt" "as.matrix(X[, index2])qsec" "as.matrix(X[, index2])gear"
  .. ..- attr(*, "assign")= int [1:4] 0 1 1 1
  ..$ qraux: num [1:4] 1.18 1.05 1.08 1.15
  ..$ pivot: int [1:4] 1 2 3 4
  ..$ tol  : num 1e-07
  ..$ rank : int 4
  ..- attr(*, "class")= chr "qr"
 $ df.residual  : int 28
 $ xlevels      : Named list()
 $ call         : language lm(formula = y ~ as.matrix(X[, index2]))
 $ terms        :Classes 'terms', 'formula' length 3 y ~ as.matrix(X[, index2])
  .. ..- attr(*, "variables")= language list(y, as.matrix(X[, index2]))
  .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:2] "y" "as.matrix(X[, index2])"
  .. .. .. ..$ : chr "as.matrix(X[, index2])"
  .. ..- attr(*, "term.labels")= chr "as.matrix(X[, index2])"
  .. ..- attr(*, "order")= int 1
  .. ..- attr(*, "intercept")= int 1
  .. ..- attr(*, "response")= int 1
  .. ..- attr(*, ".Environment")=<environment: 0x2636140> 
  .. ..- attr(*, "predvars")= language list(y, as.matrix(X[, index2]))
  .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "nmatrix.3"
  .. .. ..- attr(*, "names")= chr [1:2] "y" "as.matrix(X[, index2])"
 $ model        :'data.frame':	32 obs. of  2 variables:
  ..$ y                     : num [1:32] 21 21 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 ...
  ..$ as.matrix(X[, index2]): num [1:32, 1:3] 2.62 2.88 2.32 3.21 3.44 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:32] "Mazda RX4" "Mazda RX4 Wag" "Datsun 710" "Hornet 4 Drive" ...
  .. .. ..$ : chr [1:3] "wt" "qsec" "gear"
  ..- attr(*, "terms")=Classes 'terms', 'formula' length 3 y ~ as.matrix(X[, index2])
  .. .. ..- attr(*, "variables")= language list(y, as.matrix(X[, index2]))
  .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2] "y" "as.matrix(X[, index2])"
  .. .. .. .. ..$ : chr "as.matrix(X[, index2])"
  .. .. ..- attr(*, "term.labels")= chr "as.matrix(X[, index2])"
  .. .. ..- attr(*, "order")= int 1
  .. .. ..- attr(*, "intercept")= int 1
  .. .. ..- attr(*, "response")= int 1
  .. .. ..- attr(*, ".Environment")=<environment: 0x2636140> 
  .. .. ..- attr(*, "predvars")= language list(y, as.matrix(X[, index2]))
  .. .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "nmatrix.3"
  .. .. .. ..- attr(*, "names")= chr [1:2] "y" "as.matrix(X[, index2])"
 - attr(*, "class")= chr "lm"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(result)
## maybe str(result) ; plot(result) ...
}
\keyword{datasets}
