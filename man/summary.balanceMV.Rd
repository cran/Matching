\name{summary.balanceMV}
\alias{summary.balanceMV}
\title{Summarizing output from balanceMV}
\description{
  \code{\link{summary}} method for class \code{\link{balanceMV}}
}
\usage{
  \method{summary}{balanceMV}(object, ..., digits=5)
}
\arguments{
  \item{object}{An object of class "\code{balanceMV}", usually, a
    result of a call to \code{\link{balanceMV}}.}
  \item{digits}{The number of significant digits that should be
    displayed.}
  \item{...}{Other options for the generic summary function.}  
}
\author{
  Jasjeet S. Sekhon, Harvard University, \email{jasjeet\_sekhon@harvard.edu},
  \url{http://jsekhon.fas.harvard.edu/}
}
\seealso{
  Also see \code{\link{balanceMV}}, \code{\link{balanceUV}},
  \code{\link{summary.balanceUV}}, \code{\link{Match}},
  \code{\link{MatchBalance}}, \code{\link{GerberGreenImai}},
  \code{\link{lalonde}}
}
\examples{
data(lalonde)

#propensity score matching
glm1  <- glm(treat~age + I(age^2) + educ + I(educ^2) + black +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75, family=binomial, data=lalonde)

X  <- cbind(lalonde$re74,lalonde$re75,lalonde$age,lalonde$hisp,lalonde$black)
Y  <- lalonde$re78
Tr  <- lalonde$treat

rr  <- Match(Y=Y,Tr=Tr,X=X,M=1)

#multivariate test for balance
# 'nboots' and 'nmc' are set to small values in the interest of speed.
# Please increase to at least 500 each for publication quality p-values.  
ks  <- balanceMV(treat~ age + I(age^2) + educ + I(educ^2) + black +
               hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
               u74 + u75, data=lalonde, match.out=rr, nboots=10, nmc=10)

#let's summarize the output
summary(ks)
}
\keyword{htest}
