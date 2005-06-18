\name{GerberGreenImai}
\alias{GerberGreenImai}
\docType{data}
\title{Gerber and Green Dataset used by Imai}
\description{

  This is the dataset used by Imai (2005) to replicate and evaluate the
  field experiment done by Gerber and Green (2000). The accompanying
  demo replicates Imai's propensity score model which is then used to
  estimate the causal effect of get-out-the-vote telephone calls on
  turnout. 
  
}
\usage{data(GerberGreenImai)}
\format{
  A data frame with 10829 observations on the following 26 variables.
  \describe{
    \item{PERSONS }{tba}
    \item{WARD    }{tba}
    \item{QUESTION}{tba}
    \item{MAILGRP }{tba}
    \item{PHONEGRP}{tba}
    \item{PERSNGRP}{tba}
    \item{APPEAL  }{tba}
    \item{CONTACT }{tba}
    \item{MAILINGS}{tba}
    \item{AGE     }{tba}
    \item{MAJORPTY}{tba}
    \item{VOTE96.0}{tba}
    \item{VOTE96.1}{tba}
    \item{MAILCALL}{tba}
    \item{VOTED98 }{tba}
    \item{PHNSCRPT}{tba}
    \item{DIS.MC  }{tba}
    \item{DIS.PHN }{tba}
    \item{PHN.C   }{tba}
    \item{PHNTRT1 }{tba}
    \item{PHNTRT2 }{tba}
    \item{PHN.C1  }{tba}
    \item{PHN.C2  }{tba}
    \item{NEW     }{tba}
    \item{phone   }{tba}
    \item{AGE2    }{tba}
  }
}
\details{

  The demo provided, entitled \code{GerberGreenImai}, uses Imai's
  propensity score model to estimate the causal effect of
  get-out-the-vote telephone calls on turnout. The propensity score
  model fails to balance age.

}
\references{

  Gerber, Alan S. and Donald P. Green. 2000. ``The Effects of
  Canvassing, Telephone Calls, and Direct Mail on Voter Turnout: A Field
  Experiment.'' \emph{American Political Science Review} 94: 653-663.
  
  Imai, Kosuke. 2005. ``Do Get-Out-The-Vote Calls Reduce Turnout?
  The Importance of Statistical Methods for Field Experiments.'' 
  \emph{American Political Science Review} 99: 283-300.
  
} 
\seealso{ Also see \code{\link{Match}} and \code{\link{MatchBalance}},
  \code{\link{GenMatch}}, \code{\link{balanceMV}},
  \code{\link{balanceUV}}, \code{\link{ks.boot}}
  \code{\link{lalonde}}
}
\keyword{datasets}

