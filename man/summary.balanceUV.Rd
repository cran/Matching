\name{summary.balanceUV}
\alias{summary.balanceUV}
\title{Summarizing output from balanceUV}
\description{
  \code{\link{summary}} method for class \code{\link{balanceUV}}
}
\usage{
  \method{summary}{balanceUV}(object, ..., digits=5)
}
\arguments{
  \item{object}{An object of class "\code{balanceUV}", usually, a
    result of a call to \code{\link{balanceUV}}.}
  \item{digits}{The number of significant digits that should be
    displayed.}
  \item{...}{Other options for the generic summary function.}  
}
\author{
  Jasjeet S. Sekhon, Harvard University, \email{jasjeet\_sekhon@harvard.edu},
  \url{http://jsekhon.fas.harvard.edu/}
}
\seealso{
  Also see \code{\link{balanceUV}}, \code{\link{balanceMV}},
  \code{\link{summary.balanceMV}}, \code{\link{Match}},
  \code{\link{MatchBalance}}, \code{\link{GerberGreenImai}},
  \code{\link{lalonde}}
}
\examples{
data(lalonde)
attach(lalonde)

foo  <- balanceUV(re75[treat==1],re75[treat!=1])
summary(foo)

foo  <- balanceUV(hisp[treat==1],hisp[treat!=1])
summary(foo)
}
\keyword{htest}
