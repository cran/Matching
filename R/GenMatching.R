#RMatchingC.reload <- function()
#  {
#    if (is.loaded("FastMatchC", PACKAGE="Matching"))
#          {
#            dyn.unload("Matching.so");
#          }
#        dyn.load("Matching.so");
#  }
#RMatchingC.reload();

FastMatchC <- function(N, xvars, All, M, cdd, ww, Tr, Xmod, weights)
  {
    ret <- .Call("FastMatchC", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                 as.double(cdd), as.real(ww), as.real(Tr),
                 as.real(Xmod), as.real(weights),
                 PACKAGE="Matching")
    return(ret)
  }

FastMatchCcaliper <- function(N, xvars, All, M, cdd, caliperflag, ww, Tr, Xmod, weights, CaliperVec, Xorig)
  {
    ret <- .Call("FastMatchCcaliper", as.integer(N), as.integer(xvars), as.integer(All), as.integer(M),
                 as.double(cdd), as.integer(caliperflag), as.real(ww), as.real(Tr),
                 as.real(Xmod), as.real(weights), as.real(CaliperVec), as.real(Xorig),
                 PACKAGE="Matching")
    return(ret)
  }


MatchGenoudStage1 <- function(Tr=Tr, X=X, All=All, M=M, weights=weights)
  {
    N  <- nrow(X)
    xvars <- ncol(X)
    
# if SATC is to be estimated the treatment indicator is reversed    
    if (All==2)
      Tr <- 1-Tr

# check on the number of matches, to make sure the number is within the limits
# feasible given the number of observations in both groups.
    if (All==1)
      {
        M <- min(M,min(sum(Tr),sum(1-Tr)));        
      } else {
        M <- min(M,sum(1-Tr));
      }

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
    Mu.X  <- matrix(0, xvars, 1)
    Sig.X <- matrix(0, xvars, 1)

    weights.sum <- sum(weights)

    for (k in 1:xvars)
      {
        Mu.X[k,1] <- sum(X[,k]*weights)/weights.sum;
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weights)/weights.sum-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        X[,k]=eps/Sig.X[k,1]
      } #end of k loop

    ret <- list(Tr=Tr, X=X, All=All, M=M, N=N)
    return(ret)
  } #end of MatchGenoudStage1


FastMatchGenoud <- function(Tr, X, All=1, M=1, Weight=1, 
                            Weight.matrix=NULL,
                            weights=rep(1,length(Tr)),
                            tolerance=0.00001,
                            distance.tolerance=0.00001)
  {
    N  <- nrow(X)
    xvars <- ncol(X)

    if (!is.null(Weight.matrix))
      {
        if (Weight==2)
          {
            warning("User supplied 'Weight.matrix' is being used even though 'Weight' is not set equal to 3")
          }
        Weight  <- 3
      } else {
        Weight.matrix <- diag(ncol(X))
      }    
    
# if SATC is to be estimated the treatment indicator is reversed    
    if (All==2)
      Tr <- 1-Tr

# check on the number of matches, to make sure the number is within the limits
# feasible given the number of observations in both groups.
    if (All==1)
      {
        M <- min(M,min(sum(Tr),sum(1-Tr)));        
      } else {
        M <- min(M,sum(1-Tr));
      }

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
    Mu.X  <- matrix(0, xvars, 1)
    Sig.X <- matrix(0, xvars, 1)

    weights.sum <- sum(weights)
    
    for (k in 1:xvars)
      {
        Mu.X[k,1] <- sum(X[,k]*weights)/weights.sum;
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weights)/weights.sum-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        X[,k]=eps/Sig.X[k,1]
      } #end of k loop

    if (Weight==1)
      {
        Weight.matrix=diag(xvars)
      } else if (Weight==2) {
        if (min (eigen( t(X)%*%X/N, only.values=TRUE)$values) < tolerance)
          {
            Weight.matrix= solve(t(X)%*%X/N) 
          } else {
            Weight.matrix <- diag(xvars)
          }
      }    

    if ( min(eigen(Weight.matrix, only.values=TRUE)$values) < tolerance )
      Weight.matrix <- Weight.matrix + diag(xvars)*tolerance
        
    ww <- chol(Weight.matrix)

    rr <- FastMatchC(N=N, xvars=xvars, All=All, M=M,
                     cdd=distance.tolerance, ww=ww, Tr=Tr, Xmod=X,
                     weights=weights)    

    return(rr)
  } #end of 


###############################################################################
## For Caliper!
##
###############################################################################

MatchGenoudStage1caliper <- function(Tr=Tr, X=X, All=All, M=M, weights=weights,
                                     exact=exact, caliper=caliper)
  {
    Xorig  <- X;
    N  <- nrow(X)
    xvars <- ncol(X)

    if (!is.null(exact))
      {
        exact = as.vector(exact)
        nexacts = length(exact)
        if ( (nexacts > 1) & (nexacts != xvars) )
          {
            warning("length of exact != ncol(X). Ignoring exact option")
            exact <- NULL
          } else if (nexacts==1 & (xvars > 1) ){
            exact <- rep(exact, xvars)
          }
      }

    if (!is.null(caliper))
      {
        caliper = as.vector(caliper)
        ncalipers = length(caliper)
        if ( (ncalipers > 1) & (ncalipers != xvars) )
          {
            warning("length of caliper != ncol(X). Ignoring caliper option")
            caliper <- NULL
          } else if (ncalipers==1 & (xvars > 1) ){
            caliper <- rep(caliper, xvars)
          }
      }

    if (!is.null(caliper))
      {
        ecaliper <- vector(mode="numeric", length=xvars)
        sweights  <- sum(weights.orig)
        for (i in 1:xvars)
          {
            meanX  <- sum( X[,i]*weights.orig )/sweights
            sdX  <- sqrt(sum( (X[,i]-meanX)^2 )/sweights)
            ecaliper[i]  <- caliper[i]*sdX
          }
      } else {
        ecaliper <- NULL
      }

    if (!is.null(exact))
      {
        if(is.null(caliper))
          {
            max.diff <- abs(max(X)-min(X) + tolerance * 100)
            ecaliper <- matrix(max.diff, nrow=xvars, ncol=1)
          }
        
        for (i in 1:xvars)
          {
            if (exact[i])
              ecaliper[i] <- tolerance;
          }
      }        
    
# if SATC is to be estimated the treatment indicator is reversed    
    if (All==2)
      Tr <- 1-Tr

# check on the number of matches, to make sure the number is within the limits
# feasible given the number of observations in both groups.
    if (All==1)
      {
        M <- min(M,min(sum(Tr),sum(1-Tr)));        
      } else {
        M <- min(M,sum(1-Tr));
      }

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
    Mu.X  <- matrix(0, xvars, 1)
    Sig.X <- matrix(0, xvars, 1)

    weights.sum <- sum(weights)

    for (k in 1:xvars)
      {
        Mu.X[k,1] <- sum(X[,k]*weights)/weights.sum;
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weights)/weights.sum-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        X[,k]=eps/Sig.X[k,1]
      } #end of k loop

    ret <- list(Tr=Tr, X=X, All=All, M=M, N=N, ecaliper=ecaliper)
    return(ret)
  } #end of MatchGenoudStage1caliper


FastMatchGenoudCaliper <- function(Tr, X, All=1, M=1, Weight=1, 
                                   Weight.matrix=NULL,
                                   weights=rep(1,length(Tr)),
                                   exact=NULL, caliper=NULL,
                                   tolerance=0.00001,
                                   distance.tolerance=0.00001)
  {

    Xorig  <- X;
    N  <- nrow(X)
    xvars <- ncol(X)

    if (!is.null(Weight.matrix))
      {
        if (Weight==2)
          {
            warning("User supplied 'Weight.matrix' is being used even though 'Weight' is not set equal to 3")
          }
        Weight  <- 3
      } else {
        Weight.matrix <- diag(xvars)
      }

    if (!is.null(exact))
      {
        exact = as.vector(exact)
        nexacts = length(exact)
        if ( (nexacts > 1) & (nexacts != xvars) )
          {
            warning("length of exact != ncol(X). Ignoring exact option")
            exact <- NULL
          } else if (nexacts==1 & (xvars > 1) ){
            exact <- rep(exact, xvars)
          }
      }

    if (!is.null(caliper))
      {
        caliper = as.vector(caliper)
        ncalipers = length(caliper)
        if ( (ncalipers > 1) & (ncalipers != xvars) )
          {
            warning("length of caliper != ncol(X). Ignoring caliper option")
            caliper <- NULL
          } else if (ncalipers==1 & (xvars > 1) ){
            caliper <- rep(caliper, xvars)
          }
      }

    if (!is.null(caliper))
      {
        ecaliper <- vector(mode="numeric", length=xvars)
        sweights  <- sum(weights.orig)
        for (i in 1:xvars)
          {
            meanX  <- sum( X[,i]*weights.orig )/sweights
            sdX  <- sqrt(sum( (X[,i]-meanX)^2 )/sweights)
            ecaliper[i]  <- caliper[i]*sdX
          }
      } else {
        ecaliper <- NULL
      }

    if (!is.null(exact))
      {
        if(is.null(caliper))
          {
            max.diff <- abs(max(X)-min(X) + tolerance * 100)
            ecaliper <- matrix(max.diff, nrow=xvars, ncol=1)
          }
        
        for (i in 1:xvars)
          {
            if (exact[i])
              ecaliper[i] <- tolerance;
          }
      }            
    
# if SATC is to be estimated the treatment indicator is reversed    
    if (All==2)
      Tr <- 1-Tr

# check on the number of matches, to make sure the number is within the limits
# feasible given the number of observations in both groups.
    if (All==1)
      {
        M <- min(M,min(sum(Tr),sum(1-Tr)));        
      } else {
        M <- min(M,sum(1-Tr));
      }

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
    Mu.X  <- matrix(0, xvars, 1)
    Sig.X <- matrix(0, xvars, 1)

    weights.sum <- sum(weights)
    
    for (k in 1:xvars)
      {
        Mu.X[k,1] <- sum(X[,k]*weights)/weights.sum;
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weights)/weights.sum-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        X[,k]=eps/Sig.X[k,1]
      } #end of k loop

    if (Weight==1)
      {
        Weight.matrix=diag(xvars)
      } else if (Weight==2) {
        if (min (eigen( t(X)%*%X/N, only.values=TRUE)$values) < tolerance)
          {
            Weight.matrix= solve(t(X)%*%X/N) 
          } else {
            Weight.matrix <- diag(xvars)
          }
      }    

    if ( min(eigen(Weight.matrix, only.values=TRUE)$values) < tolerance )
      Weight.matrix <- Weight.matrix + diag(xvars)*tolerance
        
    ww <- chol(Weight.matrix)

    if (is.null(ecaliper))
      {
        caliperFlag  <- 0
        Xorig  <- 0
        CaliperVec  <- 0
      } else {
        caliperFlag  <- 1
        CaliperVec  <- ecaliper
      }    

    rr <- FastMatchC(N=N, xvars=xvars, All=All, M=M, cdd=distance.tolerance,
                     caliperflag=caliperFlag,                     
                     ww=ww, Tr=Tr, Xmod=X, weights=weights,
                     CaliperVec=CaliperVec, Xorig=Xorig)    

    return(rr)
  } #end of

###############################################################################
## GenMatchCaliper
##
###############################################################################


GenMatchCaliper <- function(Tr, X, BalanceMatrix=X, estimand="ATT", M=1,
                            weights=rep(1,length(Tr)),
                            pop.size = 50, max.generations=100,
                            wait.generations=4, hard.generation.limit=0,
                            starting.values=rep(1,ncol(X)),
                            data.type.integer=TRUE,
                            MemoryMatrix=TRUE,
                            exact=NULL, caliper=NULL,
                            nboots=0, ks=TRUE, verbose=FALSE,
                            tolerance=0.000001,
                            distance.tolerance=tolerance,
                            min.weight=0,
                            max.weight=1000,
                            Domains=NULL,
                            print.level=print.level, ...)
  {
    
    nvars <- ncol(X)
    balancevars <- ncol(BalanceMatrix)
    weighting <- 3;

    if (is.null(Domains))
      {
        Domains <- matrix(min.weight, nrow=nvars, ncol=2)
        Domains[,2] <- max.weight
      }
    
    isunix  <- Sys.getenv("OSTYPE")=="linux" | Sys.getenv("OSTYPE")=="darwin"
    if (print.level < 3 & isunix)
      {
        project.path="/dev/null"
      } else {
        project.path=paste(tempdir(),"/genoud.pro",sep="")

        #work around for rgenoud bug
        if (print.level==3)
          print.level <- 2
      }
    
    # create All
    if (estimand=="ATT")
      {
        All  <- 0
      } else if(estimand=="ATE") {
        All  <- 1
      } else if(estimand=="ATC") {
        All  <- 2
      } else {
        All  <- 0
        warning("estimand is not valid.  Defaulting to 'ATT'")
      }

    #stage 1 Match, only needs to be called once
    s1 <- MatchGenoudStage1caliper(Tr=Tr, X=X, All=All, M=M, weights=weights,
                                   exact=exact, caliper=caliper);
    s1.Tr <- s1$Tr
    s1.X <- s1$X
    s1.All <- s1$All
    s1.M <- s1$M
    s1.N <- s1$N
    s1.ecaliper  <- s1$ecaliper
    
    if (is.null(ecaliper))
      {
        caliperFlag  <- 0
        Xorig  <- 0
        CaliperVec  <- 0
      } else {
        caliperFlag  <- 1
        Xorig  <- X
        CaliperVec  <- s1$ecaliper
      }
    rm(s1)
    
    diag.nvars.tol <- diag(nvars)*tolerance
    
    genoudfunc  <- function(x)
      {
        wmatrix <- diag(x)
        if ( min(eigen(wmatrix, only.values=TRUE)$values) < tolerance )
          wmatrix <- wmatrix + diag.nvars.tol
        
        ww <- chol(wmatrix)

        rr <- FastMatchC(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                         cdd=distance.tolerance,
                         caliperflag=caliperFlag,
                         ww=ww, Tr=s1.Tr, Xmod=s1.X, weights=weights,
                         CaliperVec=CaliperVec, Xorig=Xorig)
        
        a <- GenBalance(rr=rr, X=BalanceMatrix, nvars=balancevars, nboots=nboots,
                        ks=ks, verbose=verbose)
        
        return(a)
      }    
    
    rr <- genoud(genoudfunc, nvars=nvars, starting.values=starting.values,
                 pop.size=pop.size, max.generations=max.generations,
                 wait.generations=wait.generations, hard.generation.limit=hard.generation.limit,
                 Domains=Domains,
                 MemoryMatrix=MemoryMatrix,
                 max=TRUE, gradient.check=FALSE, data.type.int=data.type.integer,
                 hessian=FALSE,
                 BFGS=FALSE, project.path=project.path,
                 print.level=print.level, ...)

    wmatrix <- diag(rr$par)
    if ( min(eigen(wmatrix, only.values=TRUE)$values) < tolerance )
      wmatrix <- wmatrix + diag.nvars.tol
        
    ww <- chol(wmatrix)
    
    mout <- FastMatchCaliper(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                             cdd=distance.tolerance,
                             caliperflag=caliperFlag,
                             ww=ww, Tr=s1.Tr, Xmod=s1.X, weights=weights,
                             CaliperVec=CaliperVec, Xorig=Xorig)

    rr2 <- list(value=rr$value, par=rr$par, Weight.matrix=wmatrix, matches=mout, ecaliper=CaliperVec)

    class(rr2) <- "GenMatch"    
    return(rr2)
  } #end of GenMatchCaliper




###############################################################################
## GenMatch
##
###############################################################################

GenMatch <- function(Tr, X, BalanceMatrix=X, estimand="ATT", M=1,
                     weights=rep(1,length(Tr)),
                     pop.size = 50, max.generations=100,
                     wait.generations=4, hard.generation.limit=FALSE,
                     starting.values=rep(1,ncol(X)),
                     data.type.integer=TRUE,
                     MemoryMatrix=TRUE,
                     exact=NULL, caliper=NULL, 
                     nboots=0, ks=TRUE, verbose=FALSE,
                     tolerance=0.000001,
                     distance.tolerance=tolerance,
                     min.weight=0,
                     max.weight=1000,
                     Domains=NULL,
                     print.level=2, ...)
  {

    Tr <- as.matrix(Tr)
    X  <- as.matrix(X)
    BalanceMatrix  <- as.matrix(BalanceMatrix)

    if(!is.null(caliper) | !is.null(exact))
      {
        rr <- GenMatchCaliper(Tr=Tr, X=X, BalanceMatrix=BalanceMatrix, estimand=estimand, M=M,
                              weights=weights,
                              caliper=caliper, exact=exact,
                              nboots=nboots, ks=ks, verbose=verbose,
                              pop.size=pop.size, max.generations=max.generations,
                              wait.generations=wait.generations,
                              hard.generation.limit=hard.generation.limit,
                              starting.values=starting.values,
                              data.type.integer=data.type.integer,
                              MemoryMatrix=MemoryMatrix,
                              tolerance=tolerance,
                              distance.tolerance=distance.tolerance,
                              max.weight=max.weight,
                              print.level=print.level, ...)
        return(rr)
      } #!is.null(caliper) | !is.null(exact)

    isunix  <- Sys.getenv("OSTYPE")=="linux" | Sys.getenv("OSTYPE")=="darwin"
    if (print.level < 3 & isunix)
      {
        project.path="/dev/null"
      } else {
        project.path=paste(tempdir(),"/genoud.pro",sep="")

        #work around for rgenoud bug
        if (print.level==3)
          print.level <- 2
      }
    
    nvars <- ncol(X)
    balancevars <- ncol(BalanceMatrix)
    weighting <- 3;

    if (is.null(Domains))
      {
        Domains <- matrix(min.weight, nrow=nvars, ncol=2)
        Domains[,2] <- max.weight
      }

    # create All
    if (estimand=="ATT")
      {
        All  <- 0
      } else if(estimand=="ATE") {
        All  <- 1
      } else if(estimand=="ATC") {
        All  <- 2
      } else {
        All  <- 0
        warning("estimand is not valid.  Defaulting to 'ATT'")
      }

    #stage 1 Match, only needs to be called once
    s1 <- MatchGenoudStage1(Tr=Tr, X=X, All=All, M=M, weights=weights);
    s1.Tr <- s1$Tr
    s1.X <- s1$X
    s1.All <- s1$All
    s1.M <- s1$M
    s1.N <- s1$N
    rm(s1)

    diag.nvars.tol <- diag(nvars)*tolerance

    genoudfunc  <- function(x)
      {
        wmatrix <- diag(x)
        if ( min(eigen(wmatrix, only.values=TRUE)$values) < tolerance )
          wmatrix <- wmatrix + diag.nvars.tol
        
        ww <- chol(wmatrix)

        rr <- FastMatchC(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                         cdd=distance.tolerance, ww=ww, Tr=s1.Tr, Xmod=s1.X,
                         weights=weights)
        
        a <- GenBalance(rr=rr, X=BalanceMatrix, nvars=balancevars, nboots=nboots,
                        ks=ks, verbose=verbose)
        
        return(a)
      }    
    
    rr <- genoud(genoudfunc, nvars=nvars, starting.values=starting.values,
                 pop.size=pop.size, max.generations=max.generations,
                 wait.generations=wait.generations, hard.generation.limit=hard.generation.limit,
                 Domains=Domains,
                 MemoryMatrix=MemoryMatrix,
                 max=TRUE, gradient.check=FALSE, data.type.int=data.type.integer,
                 hessian=FALSE,
                 BFGS=FALSE, project.path=project.path, print.level=print.level, ...)

    wmatrix <- diag(rr$par)
    if ( min(eigen(wmatrix, only.values=TRUE)$values) < tolerance )
      wmatrix <- wmatrix + diag.nvars.tol
        
    ww <- chol(wmatrix)
    
    mout <- FastMatchC(N=s1.N, xvars=nvars, All=s1.All, M=s1.M,
                       cdd=distance.tolerance, ww=ww, Tr=s1.Tr, Xmod=s1.X,
                       weights=weights)


    rr2 <- list(value=rr$value, par=rr$par, Weight.matrix=wmatrix, matches=mout, ecaliper=NULL)
    class(rr2) <- "GenMatch"
    return(rr2)
  } #end of GenMatch
