\name{summary.Match}
\alias{summary.Match}
\title{Summarizing output from Match}
\description{
  \code{\link{summary}} method for class \code{\link{Match}}
}
\usage{
  \method{summary}{Match}(object, ... , full=FALSE, digits=5)
}
\arguments{
  \item{object}{An object of class "\code{Match}", usually, a
    result of a call to \code{\link{Match}}.}
  \item{full}{A flag for whether the unadjusted estimates and naive
    standard errors should also be summarized.}
  \item{digits}{The number of significant digits that should be
    displayed.}
  \item{...}{Other options for the generic summary function.}
}
\author{
  Jasjeet S. Sekhon, Harvard University, \email{jasjeet\_sekhon@harvard.edu},
  \url{http://jsekhon.fas.harvard.edu/}
}
\seealso{ Also see \code{\link{Match}}, \code{\link{MatchBalance}}, \code{\link{balanceMV}},
  \code{\link{balanceUV}}, \code{\link{GerberGreenImai}}, \code{\link{lalonde}} }
\examples{
#
# Replication of Dehejia and Wahba psid3 model
#
# Dehejia, Rajeev and Sadek Wahba. 1999.``Causal Effects in Non-Experimental Studies: Re-Evaluating the
# Evaluation of Training Programs.''Journal of the American Statistical Association 94 (448): 1053-1062.
#
data(lalonde)

#
# Estimate the propensity model
#
glm1  <- glm(treat~age + I(age^2) + educ + I(educ^2) + black +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75, family=binomial, data=lalonde)


#
#save data objects
#
X  <- glm1$fitted
Y  <- lalonde$re78
Tr  <- lalonde$treat

#
# one-to-one matching with replacement (the "M=1" option).
# Estimating the treatment effect on the treated (the "estimand" option which defaults to 0).
#
rr  <- Match(Y=Y,Tr=Tr,X=X,M=1);

#Let's summarize the output
summary(rr)
}
\keyword{htest}
