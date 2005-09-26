\name{balanceMV}
\alias{balanceMV}
\title{Multivariate Tests for Balance}
\description{
  This function conducts the Kolmogorov-Smirnov tests for balance. The
  function does no matching.  Matching is performed by the
  \code{\link{Match}} function.  \code{balanceMV} is
  used to determine if \code{\link{Match}} was successful in achieving
  balance.  Multivariate balance is determined by the use of a model.
  This approach can be used regardless of the algorithm used to do the
  original matching.  For example, even if \code{\link{Match}} was
  told NOT to use a propensity score, \code{balanceMV} can be used to
  test balance.  Output can be summarized by using the
  \code{\link{summary.balanceMV}} function. Generally, users should call
  \code{\link{MatchBalance}} and not this function directly.
}
\usage{
balanceMV(formul, data = NULL, match.out=NULL, maxit = 1000,
          weights = rep(1,nrow(data)), nboots=100, nmc=nboots, verbose=0, ...)
}
\arguments{
  \item{formul}{ A formula denoting the model for which balance should
    be determined.  This model will be estimated by a binary logistic
    estimator. The dependent variable in the formula is usually the
    treatment indicator. Univariate balance tests will be conducted for
    each of the regressors included in this model. And the multivariate
    tests will be conducted on the predicted probabilities of treatment
    for both treated and control based on this formula. The predicted
    probability densities for both treated and control should be
    indistinguishable if balance has been achieved. Note that the model
    defined by this formula is estimated separately for the matched and
    unmatched datasets.}
  \item{data}{ A data frame which contains all of the variables in the
    formula.  If a data frame is not provided, the variables are obtained
    via lexical scoping.}
  \item{match.out}{ The output object from the \code{\link{Match}}
    function.  If this output is included, \code{balanceMV} will provide
    balance tests for both before and after matching. Otherwise
    balance tests will only be conducted for the raw data.}
  \item{maxit}{ The maximum number of iterations for the glm logistic
    procedure. }
  \item{weights}{ A vector of observation weights. }
  \item{nboots}{The number of bootstrap samples to be run.  If zero, no
    bootstraps are done.  Bootstrapping is highly recommended because
    the Kolmogorov-Smirnov test only provides correct coverage when
    bootstrapped due to the existence of nuisance parameters.  At least
    500 \code{nboots} (preferably 1000) are recommended for publication
    quality p-values.  Also see the \code{nmc} option.}
  \item{nmc}{The number of Monte Carlo simulations to be conducted for
    each Kolmogorov-Smirnov test calculated.  Monte Carlo simulations
    are highly recommended because the usual Kolmogorov-Smirnov test is not
    consistent when the densities being compared contain point masses.  At least
    500 \code{nmc} (preferably 1000) are recommended for publication
    quality p-values. Also see the \code{nboots} option.}  
  \item{verbose}{The amount of printing to be done.  If zero, there is
    no printing.  If two, details are printing such as the number of the
    bootstrap sample currently being estimated.}
  \item{...}{Further arguments passed to \code{\link{glm.control}}.}  
}
\value{
  \item{pval.kboot.unmatched}{The bootstrap p-value of the Kolmogorov-Smirnov test
    for the hypothesis that the probability densities for both the
    treated and control groups are the same (unmatched data).}
  \item{pval.kboot.matched}{The bootstrap p-value of the Kolmogorov-Smirnov test
    for the hypothesis that the probability densities for both the
    treated and control groups are the same (matched data).}
  \item{logit.unmatched}{Return object from glm estimating the model
    defined by the user in \code{formul} on the \emph{unmatched} data.}
  \item{logit.matched}{Return object from glm estimating the model
    defined by the user in \code{formul} on the \emph{matched} data.}
  \item{ks.unmatched}{Return object from \code{\link{ks.test}} before matching.}
  \item{ks.matched}{Return object from \code{\link{ks.test}} after
    matching.}
  \item{data}{A list containing five variables: a treatment indicator
    for the unmatched (\code{Tr}) and the matched data
    (\code{Tr.matched}), weights for the matched data used for the glm
    models, index.treated.indata, and index.control.indata.
    The two indexes are defined in \code{\link{Match}}.  The user should
    never need to examine this object.}
}
\details{
  The function can be used to determine if matching was successful in
  achieving balance.  Two multivariate tests are provided.  The
  Kolmogorov-Smirnov (KS) and Chi-Square null deviance tests.  The KS
  test is to be preferred over the Chi-Square test because the
  Chi-Square test is not testing the relevant hypothesis.  The null
  hypothesis for the KS test is of equal balance in the estimated
  probabilities between treated and control.  The null hypothesis for
  the Chi-Square test, however, is of all of the parameters being
  insignificant; a comparison of residual versus null deviance.  If the
  covariates being considered are discrete, this KS test is
  asymptotically nonparametric as long as the logit model does not produce zero
  parameter estimates.  The bootstrap-Monte Carlo version of the KS
  test is highly recommended because the usual KS test is not consistent
  when there are point masses in the distributions being compared, and
  the bootstrap is needed because parameters are being estimated in the
  logit model.
}
\references{
  Sekhon, Jasjeet S. 2004.  ``Balance Tests for Matching
  Estimators.''  Working
  Paper. \url{http://sekhon.polisci.berkeley.edu/matching}

  Diamond, Alexis and Jasjeet S. Sekhon. 2005. ``Genetic Matching for
  Estimating Causal Effects: A General Multivariate Matching Method for
  Achieving Balance in Observational Studies.''  Working Paper.
  \url{http://sekhon.polisci.berkeley.edu/papers/GenMatch.pdf} 

  Abadie, Alberto.  2002.  ``Bootstrap Tests for Distributional Treatment
  Effects in Instrumental Variable Models.'' \emph{Journal of the
  American Statistical Association}, 97:457 (March) 284-292.
    
  Hall, Peter. 1992. \emph{The Bootstrap and Edgeworth Expansion}.  New
  York: Springer-Verlag.

  Wilcox, Rand R. 1997. \emph{Introduction to Robust Estimation}.  San
  Diego, CA: Academic Press.

  William J. Conover (1971), \emph{Practical nonparametric statistics}.
  New York: John Wiley & Sons. Pages 295-301 (one-sample
  "Kolmogorov" test), 309-314 (two-sample "Smirnov" test).

  Shao, Jun and Dongsheng Tu. 1995. \emph{The Jackknife and Bootstrap}.
  New York: Springer-Verlag.
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.polisci.berkeley.edu/}.
}
\seealso{ Also see \code{\link{summary.balanceMV}},
  \code{\link{MatchBalance}}, \code{\link{balanceUV}},
  \code{\link{ks.boot}}, \code{\link{Match}}, \code{\link{GenMatch}},
  \code{\link{GerberGreenImai}}, \code{\link{lalonde}} }
\examples{
data(lalonde)

#
#direct matching on some variable
#
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
summary(ks)
}
\keyword{htest}


%  LocalWords:  balanceMV Kolmogorov MatchBalance formul maxit nrow nboots nmc
%  LocalWords:  Univariate regressors datasets glm logit emph covariates url UC
%  LocalWords:  Abadie Edgeworth Verlag Conover Shao Dongsheng seealso GenMatch
%  LocalWords:  balanceUV GerberGreenImai lalonde cbind hisp rr nodegr htest
