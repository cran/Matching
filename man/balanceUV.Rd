\name{balanceUV}
\alias{balanceUV}
\title{Univariate Balance Tests}
\description{
  This function provides a number of univariate balance tests.
  Generally, users should call \code{\link{MatchBalance}} and not this function
  directly.
}
\usage{
balanceUV(Tr, Co, weights = rep(1, length(Co)), exact = FALSE, ks=FALSE,
          nboots = 1000)
}
\arguments{
  \item{Tr}{A vector containing the treatment observations.}
  \item{Co}{A vector containing the control observations.}
  \item{weights}{A vector containing the observation specific
    weights. Only use this option when the treatment and control
    observations are paired (as they are after matching).}
  \item{exact}{A logical flag indicating if the exact Wilcoxon test
    should be used instead of the test with a correction.  See
    \code{\link{wilcox.test}} for details.}
  \item{ks}{ A logical flag for if the univariate bootstrap
    Kolmogorov-Smirnov (KS) test should be calculated.  If the ks option
    is set to true, the univariate KS test is calculated for all
    non-dichotomous variables.  The bootstrap KS test is consistent even
    for non-continuous variables.  See \code{\link{ks.boot}} for more
    details.}
  \item{nboots}{The number of bootstrap samples to be run for the
    \code{ks} test.  If zero, no bootstraps are done.  Bootstrapping is
    highly recommended because the bootstrapped Kolmogorov-Smirnov test
    only provides correct coverage even for non-continuous covariates. At
    least 500 \code{nboots} (preferably 1000) are
    recommended for publication quality p-values.}
}
\value{
  \item{sdiff}{This is the standardized difference which is the
    mean difference between treatment and control divided by the square
    root of the average variance (of treatment and control).  The
    entire quantity is then multiplied by a 100.}
  \item{mean.Tr}{The mean of the treatment group.}
  \item{mean.Co}{The mean of the control group.}
  \item{var.Tr}{The variance of the treatment group.}
  \item{var.Co}{The variance of the control group.}
  \item{p.value}{The p-value from the two-sided weighted \code{\link{t.test}}.}
  \item{var.ratio}{var.Tr/var.Co.}
  \item{pdiscordant}{The proportion of observations which are
    discordant.}
  \item{ks}{The object returned by \code{\link{ks.boot}}.}
  \item{wc}{The object returned by \code{\link{wilcox.test}}.}
  \item{tt}{The object returned by two-sided weighted \code{\link{t.test}}.}
}
\author{
  Jasjeet S. Sekhon, Harvard University, \email{jasjeet\_sekhon@harvard.edu},
  \url{http://jsekhon.fas.harvard.edu/}
}
\references{
  Hollander, Myles and Douglas A. Wolfe. 1973. \emph{Nonparametric
    statistical inference}. New York: John Wiley & Sons.

  Sekhon, Jasjeet S. 2004.  ``Balance Tests for Matching Estimators.''  Working Paper.  
  }
  \seealso{ Also see \code{\link{summary.balanceUV}},
    \code{\link{ks.boot}}, \code{\link{Match}}, \code{\link{MatchBalance}},
    \code{\link{balanceMV}}, \code{\link{GerberGreenImai}}, \code{\link{lalonde}} }
\examples{
data(lalonde)
attach(lalonde)

foo  <- balanceUV(re75[treat==1],re75[treat!=1])
summary(foo)

}
\keyword{univar}
