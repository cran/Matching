\name{GenMatch}
\alias{GenMatch}
\title{Genetic Matching}
\description{

  This function finds optimal balance using multivariate matching where
  a genetic search algorithm determines the weight each covariate is
  given.  This function finds the optimal weight each variable should be
  given by \code{\link{Match}} so as to achieve balance.  Balance is
  determined by a variety of univariate test, mainly paired t-tests for
  dichotomous variables and univariate Kolmogorov-Smirnov (KS) test for
  multinomial and continuous variables.  The loss criterion defining
  optimal balance is determined by the \code{loss} option.  The object
  returned by \code{GenMatch} can be supplied as the
  \code{Weight.matrix} option of the \code{\link{Match}} function to
  obtain estimates.  GenMatch, via the \code{cluster} option, supports
  the use of multiple computers, CPUs or cores to perform parallel
  computations.
}
\usage{
GenMatch(Tr, X, BalanceMatrix=X, estimand="ATT", M=1,
         weights=NULL,
         pop.size = 50, max.generations=100,
         wait.generations=4, hard.generation.limit=FALSE,
         starting.values=rep(1,ncol(X)),
         fit.func="pvals",
         data.type.integer=TRUE,
         MemoryMatrix=TRUE,
         exact=NULL, caliper=NULL, 
         nboots=0, ks=TRUE, verbose=FALSE,
         tolerance = 1e-05,
         distance.tolerance=tolerance,
         min.weight=0, max.weight=1000,
         Domains=NULL, print.level=2,
         project.path=NULL,
         paired=TRUE, loss=1,
         restrict=NULL,
         cluster=FALSE, balance=TRUE, ...)
}
\arguments{
  \item{Tr}{ A vector indicating the observations which are
    in the treatment regime and those which are not.  This can either be a
    logical vector or a real vector where 0 denotes control and 1 denotes
    treatment.}
  \item{X}{ A matrix containing the variables we wish to match on.
    This matrix may contain the actual observed covariates or the
    propensity score or a combination of both.}
  \item{BalanceMatrix}{ A matrix containing the variables we wish
    achieve balance on.  This is by default equal to \code{X}, but it can
    in principle be a matrix which contains more or less variables than
    \code{X} or variables which are transformed in various ways.  See
    the examples.}
  \item{estimand}{ A character string for the estimand.  The default
    estimand is "ATT", the sample average treatment effect for the
    treated. "ATE" is the sample average treatment effect (for all), and
    "ATC" is the sample average treatment effect for the controls.}
  \item{M}{ A scalar for the number of matches which should be
    found (with replacement). The default is one-to-one matching. }
  \item{weights}{ A vector the same length as \code{Y} which
    provides observations specific weights.  If none are provides, equal
    weights of 1 for each observations are assumed.}
  \item{pop.size}{Population Size.  This is the number of individuals
    \code{\link[rgenoud]{genoud}} uses to solve the optimization problem.  See
    \code{\link[rgenoud]{genoud}} for more details.}
  \item{max.generations}{ Maximum Generations.  This is the maximum
    number of generations that \code{\link[rgenoud]{genoud}} will run when attempting to
    optimize a function.  This is a \emph{soft} limit.  The maximum
    generation limit will be binding for \code{\link[rgenoud]{genoud}} only if
    \code{hard.generation.limit} has been set equal to \emph{TRUE}.  If it
    has not been set equal to \emph{TRUE}, \code{wait.generations}
    controls when \code{\link[rgenoud]{genoud}} stops. See \code{\link[rgenoud]{genoud}} for more details.}
  \item{wait.generations}{ If there is no improvement in the objective
    function in this number of generations, \code{\link[rgenoud]{genoud}} will think that it has
    found the optimum.  The other variables controlling termination are
    \code{max.generations} and \code{hard.generation.limit}.}
  \item{hard.generation.limit}{
    This logical variable determines if the \code{max.generations}
    variable is a binding constraint for \code{\link[rgenoud]{genoud}}.  If
    \code{hard.generation.limit} is \emph{FALSE}, then \code{\link[rgenoud]{genoud}} may exceed
    the \code{max.generations} count if the objective function has
    improved within a given number of generations (determined by
    \code{wait.generations}).}
  \item{starting.values}{
    This vector equal to the number of variables in \code{X}.  This
    vector contains the starting weights each of the variables is
    given. The \code{starting.values} vector is a way for the user
    to insert \emph{one} individual into the starting population.
    \code{\link[rgenoud]{genoud}} will randomly create the other individuals.  These values
    correspond to the diagonal of the \code{Weight.matrix} as described
    in detail in the \code{\link{Match}} function.}
  \item{fit.func}{The balance metric GenMatch should optimize.  The user
    may choose from the following or provide one:\cr
  \code{pvals}: maximize the p.values uses from a variety of hypothesis
  tests.\cr
  \code{qqmean.mean}: calculate the mean standardized difference in the QQ
  plot for each variable.  Minimize the mean of these differences across
  variables.\cr
  \code{qqmean.max}:  calculate the mean standardized difference in the QQ
  plot for each variable.  Minimize the maximum of these differences across
  variables.  This uses lexical optimization.  See the \code{lexical}
  option.\cr
  \code{qqmedian.mean}: calculate the median standardized difference in the QQ
  plot for each variable.  Minimize the median of these differences across
  variables.\cr
  \code{qqmedian.max}:  calculate the median standardized difference in the QQ
  plot for each variable.  Minimize the maximum of these differences across
  variables.  This uses lexical optimization.  See the \code{lexical}
  option.\cr
  \code{qqmax.mean}: calculate the maximum standardized difference in the QQ
  plot for each variable.  Minimize the mean of these differences across
  variables.\cr
  \code{qqmax.max}:  calculate the maximum standardized difference in the QQ
  plot for each variable.  Minimize the maximum of these differences across
  variables.  This uses lexical optimization.  See the \code{lexical}
  option.\cr
  The user may provide their own \code{fit.func}.  This function needs
  to return a fit function and assume that this value will be
  minimized.  This function should expect two arguments.  The first
  being the \code{matches} object returned by \code{GenMatch}---see
  below.  And the second being a matrix which contains the variables to
  be balanced---i.e., the \code{BalanceMatrix} the user provided to
  \code{GenMatch}.}
  \item{data.type.integer}{
      By default only integer weights are considered.  If this option is
      set to \code{false}, search will be done over floating point
      weights.  This is usually an unnecessary degree of precision.}
  \item{MemoryMatrix}{
    This variable controls if \code{\link[rgenoud]{genoud}} sets up a memory matrix.  Such a
    matrix ensures that \code{\link[rgenoud]{genoud}} will request the fitness evaluation
    of a given set of parameters only once. The variable may be
    \emph{TRUE} or \emph{FALSE}.  If it is \emph{FALSE}, \code{\link[rgenoud]{genoud}}
    will be aggressive in
    conserving memory.  The most significant negative implication of
    this variable being set to \emph{FALSE} is that \code{\link[rgenoud]{genoud}} will no
    longer maintain a memory
    matrix of all evaluated individuals.  Therefore, \code{\link[rgenoud]{genoud}} may request
    evaluations which it has already previously requested.  When
    the number variables in \code{X} is large, the memory matrix
    consumes a large amount of RAM.\cr
    
    \code{\link[rgenoud]{genoud}}'s memory matrix will require \emph{significantly} less
    memory if the user sets \code{hard.generation.limit} equal
    to \emph{TRUE}.  Doing this is a good way of conserving
    memory while still making use of the memory matrix structure.}
  \item{exact}{ A logical scalar or vector for whether exact matching
    should be done.  If a logical scalar is
    provided, that logical value is applied to all covariates of
    \code{X}.  If a logical vector is provided, a logical value should
    be provided for each covariate in \code{X}. Using a logical vector
    allows the user to specify exact matching for some but not other
    variables.  When exact matches are not found, observations are
    dropped.  \code{distance.tolerance} determines what is considered to
    be an exact match. The \code{exact} option takes precedence over the
    \code{caliper} option.  Obviously, if \code{exact} matching is done
    using \emph{all} of the covariates, one should not be using
    \code{GenMatch} unless the \code{distance.tolerance} has been set
    unusually high.}
  \item{caliper}{ A scalar or vector denoting the caliper(s) which
    should be used when matching.  A caliper is the distance which is
    acceptable for any match.  Observations which are outside of the
    caliper are dropped. If a scalar caliper is provided, this caliper is
    used for all covariates in \code{X}.  If a vector of calipers is
    provided, a caliper value should be provide for each covariate in
    \code{X}. The caliper is interpreted to be in standardized units.  For example,
    \code{caliper=.25} means that all matches not equal to or within .25
    standard deviations of each covariate in \code{X} are dropped.
    The \code{ecaliper} object which is returned by \code{GenMatch} shows
    the enforced caliper on the scale of the \code{X} variables.}
  \item{nboots}{The number of bootstrap samples to be run for the
    \code{ks} test.}
  \item{ks}{ A logical flag for if the univariate bootstrap
    Kolmogorov-Smirnov (KS) test should be calculated.  If the ks option
    is set to true, the univariate KS test is calculated for all
    non-dichotomous variables.  The bootstrap KS test is consistent even
    for non-continuous variables.  If a given variable is dichotomous, a
    t-test is used even if the KS test is requested.  See
    \code{\link{ks.boot}} for more details.}
  \item{verbose}{If details should be printed of each fit evaluation
    done by the genetic algorithm.  Verbose is set to FALSE if the
    \code{cluster} option is used.}
  \item{tolerance}{ This is a scalar which is used to determine
    numerical tolerances.  This option is used by numerical routines
    such as those used to determine if a matrix is singular.}
  \item{distance.tolerance}{This is a scalar which is used to determine if distances
    between two observations are different from zero.  Values less than
    \code{distance.tolerance} are deemed to be equal to zero.  This
    option can be used to perform a type of optimal matching}
  \item{min.weight}{This is the minimum weight any variable may be
    given.}  
  \item{max.weight}{This is the maximum weight any variable may be
    given.}
  \item{Domains}{This is a \code{ncol(X)} \eqn{\times 2}{*2} matrix.
    The first column is the lower bound, and the second column is the
    upper bound for each variable over which \code{\link[rgenoud]{genoud}} will
    search for weights.  If the user does not provide this matrix, the
    bounds for each variable will be determined by the \code{min.weight}
    and \code{max.weight} options.}
  \item{print.level}{
    This option controls the level of printing.  There
    are four possible levels: 0 (minimal printing), 1 (normal), 2
    (detailed), and 3 (debug).  If level 2 is selected, \code{GenMatch} will
    print details about the population at each generation, including the
    best individual found so far. If debug
    level printing is requested, details of the \code{\link[rgenoud]{genoud}}
    population are printed in the "genoud.pro" file which is located in
    the temporary \code{R} directory returned by the \code{\link{tempdir}}
    function.  See the \code{project.path} option for more details.
    Because \code{GenMatch} runs may take a long time, it is important for the
    user to receive feedback.  Hence, print level 2 has been set as the
    default.}
  \item{project.path}{ This is the path of the
    \code{\link[rgenoud]{genoud}} project file.  By default no file is
    produced unless \code{print.level=3}.  In that case,
    \code{\link[rgenoud]{genoud}} places it's output in a file called
    "genoud.pro" located in the temporary directory provided by
    \code{\link{tempdir}}.  If a file path is provided to the
    \code{project.path} option, a file will be created regardless of the
    \code{print.level}. The behavior of the project file, however, will
    depend on the \code{print.level} chosen.  If the \code{print.level}
    variable is set to 1, then the project file is rewritten after each
    generation.  Therefore, only the currently fully completed generation
    is included in the file.  If the \code{print.level} variable is set to
    2 or higher, then each new generation is simply appended to the
    project file. No project file is generated for
    \code{print.level=0}.}
  \item{paired}{A flag for if the paired \code{\link{t.test}} should be
    used when determining balance.}
  \item{loss}{The loss function to be optimized.  The default value, \code{1},
    implies "lexical" optimization: all of the balance statistics will
    be sorted from the most discrepant to the least and weights will be
    picked which minimize the maximum discrepancy.  If multiple sets of
    weights result in the same maximum discrepancy, then the second
    largest discrepancy is examined to choose the best weights.  The
    processes continues iteratively until ties are broken.  \cr

    If the value of \code{2} is used, then only the maximum discrepancy
    is examined.  This was the default behavior prior to version 1.0.  The
    user may also pass in any function she desires.  Note that the
    option 1 corresponds to the \code{\link{sort}} function and option 2
    to the \code{\link{min}} function.  Any user specified function
    should expect a vector of balance statistics  ("p-values") and it
    should return either a vector of values (in which case "lexical"
    optimization will be done) or a scalar value (which will be
    maximized). Some possible alternative functions are
    \code{\link{mean}} or \code{\link{median}}.}
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
  \item{cluster}{This can either be an object of the  'cluster' class
    returned by one of the \code{\link[snow]{makeCluster}} commands in
    the snow package or a vector of machine names so \code{GenMatch} can
    setup the cluster automatically.  If it is the later, the vector should
    look like: \cr
    \code{c("localhost","musil","musil","deckard")}.\cr
    This vector would create a cluster with four nodes: one on the
    localhost another on "deckard" and two on the machine named
    "musil".  Two nodes on a given machine make sense if the machine has two or more
    chips/cores.  \code{GenMatch} will setup a SOCK cluster by a call to
    \code{\link[snow]{makeSOCKcluster}}.  This will require the user
    to type in her password for each node as the cluster is by default
    created via \code{ssh}.  One can add on usernames to the machine
    name if it differs from the current shell: "username@musil".  Other
    cluster types, such as PVM and MPI,
    which do not require passwords can be created by directly calling
    \code{\link[snow]{makeCluster}}, and then passing the returned
    cluster object to \code{GenMatch}. For an example of how to manually setup up
    a cluster with a direct call to \code{\link[snow]{makeCluster}} see
    \url{http://sekhon.berkeley.edu/matching/R/cluster_manual.R}.
    For an example of how to get around a firewall by ssh tunneling see:
    \url{http://sekhon.berkeley.edu/matching/R/cluster_manual_tunnel.R}.}
  \item{balance}{This logical flag controls if load balancing is done
    across the cluster.  Load balancing can result in better cluster
    utilization; however, increased communication can reduce
    performance.  This options is best used if each individual call to
    \code{\link[Matching]{Match}} takes at least several minutes to
    calculate or if the
    nodes in the cluster vary significantly in their performance. If
    cluster==FALSE, this option has no effect.}  
  \item{...}{Other options which are passed on to \code{\link[rgenoud]{genoud}}.}
}
\value{
  \item{value}{The lowest p-value of the matched dataset.}
  \item{par}{A vector of the weights given to each variable in \code{X}.}
  \item{Weight.matrix}{A matrix whose diagonal corresponds to the weight
    given to each variable in \code{X}.  This object corresponds to the
    \code{Weight.matrix} in the \code{\link{Match}} function.}
  \item{matches}{A matrix where the first column contains the row
    numbers of the treated observations in the matched dataset. The second
    column contains the row numbers of the control observations. And the
    third column contains the weight that each matched pair is given.
    These columns respectively correspond to the \code{index.treated},
    \code{index.control} and \code{weights} objects which are returned by
    \code{\link{Match}}.}
  \item{ecaliper }{The size of the enforced caliper on the scale of the
    \code{X} variables.  This object has the same length as the number of
    covariates in \code{X}.}
}
\references{
  Sekhon, Jasjeet S. 2006.  ``Alternative Balance Metrics for Bias
  Reduction in Matching Methods for Causal Inference.''  Working Paper.
  \url{http://sekhon.berkeley.edu/papers/SekhonBalanceMetrics.pdf}

  Sekhon, Jasjeet S. 2006.  ``Matching: Algorithms and Software for
  Multivariate and Propensity Score Matching with Balance Optimization
  via Genetic Search.'' \url{http://sekhon.berkeley.edu/matching/}    
  
  Diamond, Alexis and Jasjeet S. Sekhon. 2005. ``Genetic Matching for
  Estimating Causal Effects: A General Multivariate Matching Method for
  Achieving Balance in Observational Studies.''  Working Paper.
  \url{http://sekhon.berkeley.edu/papers/GenMatch.pdf}

  Sekhon, Jasjeet Singh and Walter R. Mebane, Jr. 1998.  ``Genetic
  Optimization Using Derivatives: Theory and Application to Nonlinear
  Models.'' \emph{Political Analysis}, 7: 187-210.
  \url{http://sekhon.berkeley.edu/genoud/genoud.pdf} 
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}.  
}
\seealso{ Also see \code{\link{Match}}, \code{\link{summary.Match}},
  \code{\link{MatchBalance}}, \code{\link[rgenoud]{genoud}},
  \code{\link{balanceMV}}, \code{\link{balanceUV}}, \code{\link{qqstats}},
  \code{\link{ks.boot}}, \code{\link{GerberGreenImai}}, \code{\link{lalonde}}
} 
\examples{
set.seed(38913)

data(lalonde)
attach(lalonde)

#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74);

#The covariates we want to obtain balance on
BalanceMat <- cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74,
                    I(re74*re75));

#Let's call GenMatch() to find the optimal weight to give each
#covariate in 'X' so as we have achieved balance on the covariates in
#'BalanceMat'. This is only an example so we want GenMatch to be quick
#to the population size has been set to be only 15 via the 'pop.size'
#option.  
genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATE", M=1,
                   pop.size=16, max.generations=10, wait.generations=1)

#The outcome variable
Y=re78/1000;

# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
mout <- Match(Y=Y, Tr=treat, X=X, estimand="ATE", Weight.matrix=genout)
summary(mout)

#                        
#Let's determine if balance has actually been obtained on the variables of interest
#                        
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                   re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500, ks=TRUE, mv=FALSE)

# For more examples see: http://sekhon.berkeley.edu/matching/R.
}
\keyword{nonparametric}

%  LocalWords:  GenMatch Kolmogorov multinomial CPUs BalanceMatrix estimand ATT
%  LocalWords:  ncol MemoryMatrix nboots ATC rgenoud genoud emph cr ecaliper UC
%  LocalWords:  eqn tempdir rbind makeCluster localhost musil deckard usernames
%  LocalWords:  makeSOCKcluster PVM MPI url dataset seealso MatchBalance cbind
%  LocalWords:  balanceMV balanceUV GerberGreenImai lalonde hisp nodegr genout
%  LocalWords:  BalanceMat mout mb mv pvals qqmean QQ qqmedian qqmax
