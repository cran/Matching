\name{Match}
\alias{Match}
\title{Multivariate and Propensity Score Matching Estimator for Causal Inference}
\description{

  This function performs multivariate matching.  This function is
  intended to be used in conjunction with the \code{MatchBalance}
  function which checks if the results of this function have actually
  achieved balance on a set of covariates.  If one wants to do
  propensity score matching, one should estimate the propensity model
  before calling \code{Match}, and then send \code{Match} the propensity
  scores to use. The \code{\link{GenMatch}} function can be used to
  \emph{automatically find balance} by the use of a genetic search
  algorithm which determines the optimal weight to give each covariate.
  \code{Match} provides principled standard errors when matching is done
  with covariates or a known propensity score.  Ties are handled in a
  deterministic and coherent fashion.
}
\usage{
Match(Y=NULL, Tr, X, Z = X, V = rep(1, length(Y)), estimand = "ATT", M = 1,
      BiasAdjust = FALSE, exact = NULL, caliper = NULL,
      Weight = 1, Weight.matrix = NULL, weights = rep(1, length(Y)),
      Var.calc = 0, sample = FALSE, tolerance = 1e-05,
      distance.tolerance = 1e-05, restrict=NULL,
      match.out = NULL, version="fast")
}
\arguments{
  \item{Y}{ A vector containing the outcome of interest.
    Missing values are not allowed.  An outcome vector is not required
    because the matches generated will be the same regardless of the
    outcomes.  Of course without any outcomes, no causal effect
    estimates will be produced only a matched dataset. }
  \item{Tr}{ A vector indicating the observations which are
    in the treatment regime and those which are not.  This can either be a
    logical vector or a real vector where 0 denotes control and 1 denotes
    treatment.}
  \item{X}{ A matrix containing the variables we wish to match on.
    This matrix may contain the actual observed covariates or the
    propensity score or a combination of both. All columns of this
    matrix must have positive variance or \code{Match} will return an
    error.}
  \item{Z}{ A matrix containing the covariates for which we wish to make
    bias adjustments.}
  \item{V}{ A matrix containing the covariates for which the variance
    of the causal effect may vary. Also see the \code{Var.calc} option,
    which takes precedence.}
  \item{estimand}{ A character string for the estimand.  The default
    estimand is "ATT", the sample average treatment effect for the
    treated. "ATE" is the sample average treatment effect (for all), and
    "ATC" is the sample average treatment effect for the controls.}
  \item{M}{ A scalar for the number of matches which should be
    found (with replacement). The default is one-to-one matching. }
  \item{BiasAdjust}{ A logical scalar for whether regression adjustment
    should be used. See the \code{Z} matrix.}
  \item{exact}{ A logical scalar or vector for whether exact matching
    should be done.  \emph{Variables which are to be exactly matched on
      should also be given a very large weight (e.g., 1000) via the
      \code{Weight.matrix} option}. If a logical scalar is
    provided, that logical value is applied to all covariates of
    \code{X}.  If a logical vector is provided, a logical value should
    be provided for each covariate in \code{X}. Using a logical vector
    allows the user to specify exact matching for some but not other
    variables.  When exact matches are not found, observations are
    dropped.  \code{distance.tolerance} determines what is considered to be an
    exact match. The \code{exact} option takes precedence over the
    \code{caliper} option.}
  \item{caliper}{ A scalar or vector denoting the caliper(s) which
    should be used when matching.  \emph{Variables for which a caliper is
      to be used should also be given a very large weight (e.g., 1000) via
      the \code{Weight.matrix} option}.  A caliper is the distance which is
    acceptable for any match.  Observations which are outside of the
    caliper are dropped. If a scalar caliper is provided, this caliper is
    used for all covariates in \code{X}.  If a vector of calipers is
    provided, a caliper value should be provide for each covariate in
    \code{X}. The caliper is interpreted to be in standardized units.  For
    example, \code{caliper=.25} means that all matches not equal to or
    within .25 standard deviations of each covariate in \code{X} are
    dropped.  The \code{ecaliper} object which is returned by \code{Match}
    shows the enforced caliper on the scale of the \code{X} variables.}
  \item{Weight}{ A scalar for the type of
    weighting scheme the matching algorithm should use when weighting
    each of the covariates in \code{X}.  The default value of
    1 denotes that weights are equal to the inverse of the variances. 2
    denotes the Mahalanobis distance metric, and 3 denotes
    that the user will supply a weight matrix (\code{Weight.matrix}).  Note that
    if the user supplies a \code{Weight.matrix}, \code{Weight} will be automatically
    set to be equal to 3.}
  \item{Weight.matrix}{ This matrix denotes the weights the matching
    algorithm uses when weighting each of the covariates in \code{X}---see
    the \code{Weight} option. This square matrix should have as many
    columns as the number of columns of the \code{X} matrix. This matrix
    is usually provided by a call to the \code{\link{GenMatch}} function
    which finds the optimal weight each variable should be given so as to
    achieve balance on the covariates. \cr

    For most uses, this matrix has zeros in the off-diagonal
    cells.  This matrix can be used to weight some variables more than
    others.  For
    example, if \code{X} contains three variables and we want to 
    match as best as we can on the first, the following would work well:
    \cr
    \code{> Weight.matrix <- diag(3)}\cr
    \code{> Weight.matrix[1,1] <- 1000/var(X[,1])} \cr
    \code{> Weight.matrix[2,2] <- 1/var(X[,2])} \cr
    \code{> Weight.matrix[3,3] <- 1/var(X[,3])} \cr
    This code changes the weights implied by the
    inverse of the variances by multiplying the first variable by a 1000
    so that it is highly weighted.  In order to enforce exact matching
    see the \code{exact} and \code{caliper} options.}
  \item{weights}{ A vector the same length as \code{Y} which
    provides observations specific weights.}
  \item{Var.calc}{ A scalar for the variance estimate
    that should be used.  By default \code{Var.calc=0} which means that
    homoscedasticity is assumed.  For values of  \code{Var.calc > 0},
    robust variances are calculated using \code{Var.calc} matches. }
  \item{sample}{ A logical flag for whether the population or sample
    variance is returned. }
  \item{tolerance}{ This is a scalar which is used to determine
    numerical tolerances.  This option is used by numerical routines
    such as those used to determine if a matrix is singular.}
  \item{distance.tolerance}{This is a scalar which is used to determine if distances
    between two observations are different from zero.  Values less than
    \code{distance.tolerance} are deemed to be equal to zero.  This
    option can be used to perform a type of optimal matching}
  \item{restrict}{A matrix which restricts the possible matches.  This
    matrix has one row for each restriction and three
    columns.  The first two columns contain the two observation numbers
    which are to be restricted (for example 4 and 20), and the third
    column is the restriction imposed on the observation-pair.
    Negative numbers in the third column imply that the two observations
    cannot be matched under any circumstances, and positive numbers are
    passed on as the distance between the two observations for the
    matching algorithm.  The most commonly used positive restriction is
    \code{0} which implies that the two observations will always
    be matched.  \cr

    Exclusion restriction are even more common.  For example, if we want
    to exclude the observation pair 4 and 20 and
    the pair 6 and 55 from being matched, the restrict matrix would be:
    \code{restrict=rbind(c(4,20,-1),c(6,55,-1))}}
  \item{match.out}{ The return object from a previous call to
    \code{Match}.  If this object is provided, then \code{Match} will
    use the matches found by the previous invocation of the function.
    Hence, \code{Match} will run faster.  This is
    useful when the treatment does not vary across calls to
    \code{Match} and one wants to use the same set of matches as found
    before.  This often occurs when one is trying to estimate the causal
    effect of the same treatment (\code{Tr}) on different outcomes
    (\code{Y}). Be careful when using this option to use the same
    arguments as used for the previous invocation of \code{Match} unless
    you know exactly what you are doing.}
  \item{version}{The version of the code to be used.  The "fast" C/C++
    version of the code is used unless the "old" (stable) version is
    requested.}
}
\details{
  
  This function is intended to be used in conjunction with the
  \code{MatchBalance} function which checks if the results of this
  function have actually achieved balance.  The results of this function
  can be summarized by a call to the \code{\link{summary.Match}}
  function. If one wants to do propensity score matching, one should estimate the
  propensity model before calling \code{Match}, and then place the
  fitted values in the \code{X} matrix---see the provided example.  \cr

  The \code{\link{GenMatch}} function can be used to \emph{automatically
  find balance} by the use of a genetic search algorithm which determines
  the optimal weight to give each covariate. The object returned by
  \code{\link{GenMatch}} can be supplied to the \code{Weight.matrix}
  option of \code{Match} to obtain estimates.\cr
  
  Three demos are included: \code{GerberGreenImai}, \code{DehejiaWahba},
  and \code{AbadieImbens}.  These can be run by calling the
  \code{\link{demo}} function such as by \code{demo(DehejiaWahba)}. \cr
}

\value{
  \item{est }{The estimated average causal effect.}
  \item{se }{The standard error.  This standard error is principled if
    \code{X} consists of either covariates or a known propensity score
    because it takes into account the uncertainty of the matching
    procedure.  If an estimated propensity score is used, the
    uncertainty involved in its estimation is not accounted for although the
    uncertainty of the matching procedure itself still is.}
  \item{est.noadj }{The estimated average causal effect without any
    \code{BiasAdjust}.  If \code{BiasAdjust} is not requested, this is the
    same as \code{est}.}
  \item{se.naive }{The naive standard error.  This is the standard error
    calculated on the matched data using the usual method of calculating
    the difference of means (between treated and control) weighted by the
    observation weights provided by \code{weights}.  Note that the
    standard error provided by \code{se} takes into account the uncertainty
    of the matching procedure while \code{se.naive} does not.  Neither
    \code{se} nor \code{se.naive} take into account the uncertainty of
    estimating a propensity score.  \code{se.naive} does
    not take into account any \code{BiasAdjust}.  Summary of the naive
    results can be requested by setting the \code{full=TRUE} flag when
    using the \code{\link{summary.Match}} function on the object
    returned by
    \code{Match}.}
  \item{se.cond }{The conditional standard error.  The practitioner
    should not generally use this.}
  \item{mdata }{A list which contains the matched datasets produced by
    \code{Match}.  Three datasets are included in this list: \code{Y},
    \code{Tr} and \code{X}.}
  \item{index.treated }{A vector containing the observation numbers from
    the original dataset for the treated observations in the
    matched dataset.  This index in conjunction with \code{index.control}
    can be used to recover the matched dataset produced by
    \code{Match}.  For example, the \code{X} matrix used by \code{Match}
    can be recovered by
    \code{rbind(X[index.treated,],X[index.control,])}. The user should
    generally just examine the output of \code{mdata}.}
  \item{index.control }{An index for the control observations in the
    matched data.  This index in conjunction with \code{index.treated}
    can be used to recover the matched dataset produced by
    \code{Match}.  For example, the \code{X} matrix used by \code{Match}
    can be recovered by
    \code{rbind(X[index.treated,],X[index.control,])}. The user should
    generally just examine the output of \code{mdata}.}
  \item{weights }{The weight for the matched dataset.  If all of the
    observations had a weight of 1 on input, they will have a weight of
    1 on output if each observation was only matched once.}
  \item{orig.nobs }{The original number of observations in the dataset.}
  \item{orig.wnobs }{The original number of weighted observations in the
    dataset.}
  \item{orig.treated.nobs}{The original number of treated observations (unweighted).}
  \item{nobs }{The number of observations in the matched dataset.}
  \item{wnobs }{The number of weighted observations in the matched dataset.}
  \item{caliper }{The \code{caliper} which was used.}
  \item{ecaliper }{The size of the enforced caliper on the scale of the
    \code{X} variables.  This object has the same length as the number of
    covariates in \code{X}.}
  \item{exact}{The value of the \code{exact} function argument.}
  \item{ndrops}{The number of actual observations which were dropped
    either because of caliper or exact matching.  This number is not
    reliable if observation specific weights were passed in using the \code{weights}
    option.  But \code{ndrops.matches} will always be accurate.}
  \item{ndrops.matches}{The number of matches including ties which were dropped either because of the
    caliper or exact matching. Note that since this number includes
    ties, it is not the same as \code{ndrops}.}  
}
\references{

  Abadie, Alberto and Guido Imbens. 2005.
  ``Large Sample Properties of Matching Estimators for Average
  Treatment Effects.'' \emph{Econometrica}.
  \url{http://ksghome.harvard.edu/~.aabadie.academic.ksg/sme.pdf} 

  Diamond, Alexis and Jasjeet S. Sekhon. 2005. ``Genetic Matching for
  Estimating Causal Effects: A General Multivariate Matching Method for
  Achieving Balance in Observational Studies.''  Working Paper.
  \url{http://sekhon.berkeley.edu/papers/GenMatch.pdf}

  Imbens, Guido. 2004. Matching Software for Matlab and
  Stata. \url{http://elsa.berkeley.edu/~imbens/estimators.shtml}
  
  Sekhon, Jasjeet S. 2004.  ``The Varying Role of Voter Information
  Across Democratic Societies.''  Working Paper.
  \url{http://sekhon.berkeley.edu/papers/SekhonInformation.pdf}
}
\author{Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}.
}
\seealso{ Also see \code{\link{summary.Match}},
  \code{\link{GenMatch}},
  \code{\link{MatchBalance}},
  \code{\link{balanceMV}}, \code{\link{balanceUV}}, \code{\link{ks.boot}},
  \code{\link{GerberGreenImai}}, \code{\link{lalonde}}
} 
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
summary(rr)

#
# Let's check for balance
# 'nboots' and 'nmc' are set to small values in the interest of speed.
# Please increase to at least 500 each for publication quality p-values.  
mb  <- MatchBalance(treat~age + I(age^2) + educ + I(educ^2) + black +
                    hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
                    u74 + u75, data=lalonde, match.out=rr, nboots=10, nmc=10)
}
\keyword{nonparametric}


%  LocalWords:  MatchBalance GenMatch emph estimand ATT BiasAdjust calc dataset
%  LocalWords:  ATC ecaliper cr diag homoscedasticity rbind GerberGreenImai se
%  LocalWords:  DehejiaWahba AbadieImbens noadj cond mdata datasets wnobs url
%  LocalWords:  ndrops Abadie Imbens Econometrica Matlab Stata UC seealso Wahba
%  LocalWords:  balanceMV balanceUV lalonde Dehejia psid Rajeev Sadek glm hisp
%  LocalWords:  nodegr rr nboots nmc mb
