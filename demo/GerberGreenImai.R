#
# Gerber, Alan S. and Donald P. Green. 2000. "The Effects of Canvassing, Telephone Calls, and
# Direct Mail on Voter Turnout: A Field Experiment." American Political Science Review 94: 653-663.
#
# Imai, Kosuke. forthcoming. "Do Get-Out-The-Vote Calls Reduce Turnout? The Importance of 
# Statistical Methods for Field Experiments".  http://www.princeton.edu/~kimai/research/matching.html
#

set.seed(10391)
data(GerberGreenImai)

Tr<-GerberGreenImai$PHN.C1 #treatment
Y<-GerberGreenImai$VOTED98 #outcome

Z  <- as.matrix(cbind(GerberGreenImai$WARD, GerberGreenImai$MAJORPTY,
                      GerberGreenImai$PERSONS, GerberGreenImai$VOTE96.1,
                      GerberGreenImai$NEW, GerberGreenImai$AGE, GerberGreenImai$AGE2,
                      GerberGreenImai$PERSONS*GerberGreenImai$VOTE96.1,
                      GerberGreenImai$PERSONS*GerberGreenImai$NEW))

r1  <- Match(Y=Y, Tr=Tr, Z=Z, X=Z, M=1, BiasAdj=TRUE, caliper=.1)
summary(r1, full=TRUE)
mb1  <- MatchBalance(PHN.C1 ~ PERSONS + VOTE96.1 + NEW + MAJORPTY + AGE +
                     WARD + I(PERSONS*VOTE96.1) + I(PERSONS*NEW) + AGE2, match.out=r1,
                     data=GerberGreenImai, nboots=100, nmc=500)

