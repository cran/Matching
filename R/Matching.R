# Jasjeet S. Sekhon <jasjeet_sekhon@harvard.edu>
# HTTP://jsekhon.fas.harvard.edu
# Harvard University

# Match(): function to estimate treatments using a matching estimator.
# Currently only the ability to estimate average treatment effects
# using the approach of Abadie and Imbens is implemented.  In the
# future, quantile treatment effects will be implemented along with
# the ability to use robust estimation when estimating the propensity
# score. MatchBalance(), balanceMV() and balanceUV() test for balance.

Match  <- function(Y,Tr,X,Z=X,V=rep(1,length(Y)),
                   estimand="ATT",M=1,BiasAdj=FALSE,Weight=1,Weight.matrix=NULL,
                   Var.calc=0,weights=rep(1,length(Y)), 
                   caliper=FALSE, exact=FALSE, sample=FALSE,
                   extra.output=FALSE, tolerance=0.00001)
  {
    isna  <- sum(is.na(Y)) + sum(is.na(Tr)) + sum(is.na(X)) + sum(is.na(Z))
    if (isna!=0)
      {
        stop("Match(): input includes NAs")
        return(NULL)
      }

    Y  <- as.matrix(Y)
    Tr <- as.matrix(Tr)
    X  <- as.matrix(X)
    Z  <- as.matrix(Z)
    V  <- as.matrix(V)
    weights <- as.matrix(weights)
    BiasAdj  <- as.real(BiasAdj)
    exact  <- as.real(exact)
    sample  <- as.real(sample)

    ccc  <- tolerance
    cdd  <- tolerance

    if (estimand=="ATT")
      {
        estimand  <- 0
      } else if(estimand=="ATE") {
        estimand  <- 1
      } else if(estimand=="ATC") {
        estimand  <- 2
      } else {
        estimand  <- 0
        warning("estimand is not valid.  Defaulting to 'ATT'")
      }

    if (!is.null(Weight.matrix))
      {
        if (Weight!=3)
          {
            warning("User supplied 'Weight.matrix' is being used even though 'Weight' is not set equal to 3")
          }
        Weight  <- 3
      } else {
        Weight.matrix <- dim(X)[2]
      }

    if ( (caliper < 0) )
      {
        warning("caliper < 0, ignoring option")
        caliper <- FALSE
      }

    orig.nobs  <- length(Y)
    nobs  <- orig.nobs
    orig.treated.nobs  <- sum(Tr==1)
    orig.wnobs  <- sum(weights)
    weights.orig  <- as.matrix(weights)

    ret <- Rmatch(Y=Y, Tr=Tr, X=X, Z=Z, V=V, All=estimand, M=M, BiasAdj=BiasAdj,
                  Weight=Weight, Weight.matrix=Weight.matrix, Var.calc=Var.calc,
                  weight=weights, SAMPLE=sample, exact=exact,
                  ccc=tolerance, cdd=tolerance)
    ret.use <- ret

    indx <-  cbind(ret$art.data[,1],  ret$art.data[,2],  ret$W)

    index.treated  <- indx[,1]
    index.control  <- indx[,2]
    index.treated.indata  <- index.treated
    index.control.indata  <- index.control
    weights        <- indx[,3]
    
    index.treated.nocaliper  <- index.treated
    index.control.nocaliper  <- index.control
    index.caliper  <- rep(TRUE,length(index.treated))
    index.drop  <- NULL
    
    if (caliper!=FALSE)
      {
        n  <- ncol(X)
        for (i in 1:n)
          {
            sweights  <- sum(weights.orig)
            meanX  <- sum( X[,i]*weights.orig )/sweights
            sdX  <- sqrt(sum( (X[,i]-meanX)^2 )/sweights)
            ecaliper  <- caliper*sdX
            
            index.caliper  <- (abs(X[index.treated,i]-X[index.control,i]) < ecaliper) & index.caliper
          }
#        index.drop.treated  <- as.vector(-1*(index.treated.nocaliper[!index.caliper]))
#        index.drop.control  <- as.vector(-1*(index.control.nocaliper[!index.caliper]))
#        index.drop  <- unique(c(index.drop.treated,index.drop.control))        
        index.drop  <- -1*(unique(index.treated.nocaliper[!index.caliper]))
        sum.caliper.drops  <- length(index.drop)
        
        
        if( estimand <= 1 ) {
          #ATT or ATE
          valid.matches  <- sum(Tr[index.drop]==1)
        } else {
          #ATC
          valid.matches  <- sum(Tr[index.drop]==0)
        }
          
        if (sum.caliper.drops > 0 & valid.matches < 1)
          {
            sum.caliper.drops  <- -9999
            indx  <- as.matrix(NA)
          }

        if (sum.caliper.drops > 0)
          {
            #reestimate with dropped observations
            
            cY  <- as.matrix(Y[index.drop])
            cTr <- as.matrix(Tr[index.drop])
            cweights <- as.matrix(weights.orig[index.drop])
            cX  <- as.matrix(X[index.drop,])
            cZ  <- as.matrix(Z[index.drop,])
            cV  <- as.matrix(V[index.drop,])

            ret2 <- Rmatch(Y=cY, Tr=cTr, X=cX, Z=cZ, V=cV, All=estimand, M=M, BiasAdj=BiasAdj,
                           Weight=Weight, Weight.matrix=Weight.matrix, Var.calc=Var.calc,
                           weight=cweights, SAMPLE=sample, exact=exact,
                           ccc=tolerance, cdd=tolerance)

            indx <- cbind(ret2$art.data[,1],  ret2$art.data[,2],  ret2$W)

            if(!is.na(indx[1,1]))
              {
                index.treated.indata  <- indx[,1]
                index.control.indata  <- indx[,2]
                weights        <- indx[,3]
                
                #lookup table
                vdrop  <- abs(index.drop)
                vdrop.tmp  <- vdrop
                ndrops  <- length(vdrop)
                nobs  <- orig.nobs-ndrops
                lu  <- matrix(0, nrow=nobs,ncol=1)
                d  <- 1
                count  <- 0
                i  <- 1 
                while (i <= nobs)
                  {
#                cat("i:",i,"d:",d,"vdrop.tmp[d]",vdrop.tmp[d],"\n")
                    if (d <= ndrops)
                      {
                        if (vdrop.tmp[d]==i)
                          {
                            count  <- count+2
                            lu[i]  <- count
                                        #            vdrop.tmp[d+1]  <- vdrop.tmp[d+1]-d
                            vdrop.tmp[(d+1):ndrops]  <- vdrop.tmp[(d+1):ndrops]-1
#                        cat("**IN d:",d,"\n")
#                        print(vdrop.tmp)
                            d  <- d+1
                            if (d<=ndrops)
                              {
                                if (vdrop.tmp[d]==vdrop.tmp[d-1])
                                  {
                                    i  <- i - 1
                                    count  <- count-1
#                                cat("reset i:",i,"\n")
                                  }
                              }
                          } else {
                            count  <- count+1
                            lu[i]  <- count;
                          }
                      } else {
                        count  <- count+1
                        lu[i]  <- count;        
                      }
                    i  <- i + 1
                  }
                index.treated  <- lu[index.treated.indata]
                index.control  <- lu[index.control.indata]                        
                
              }#end of !is.na(indx[1,1])
            ret.use <- ret2
          }#end of restimate for caliper portion of code
      }#end of caliper loop

    if (is.na(indx[1,1]))
      {
        if (sum.caliper.drops < 0) {
          warning("'Match' object contains no valid matches (probably because of the caliper).") 
        } else {
          warning("'Match' object contains no valid matches")
        }
      }#is.na(indx[1,1])

    if (!is.na(indx[1,1]))
      {
        if (!extra.output)
          {
            ret.use$art.data  <- NULL
            ret.use$aug.data  <- NULL
          }
        #RESET INDEX.TREATED        
        indx  <- as.matrix(cbind(index.treated,index.control))
        if (estimand==0) {
          #"ATT"
          index.treated  <- indx[,1]
          index.control  <- indx[,2]
        } else if(estimand==1) {
          #"ATE"
          tmp.index.treated  <- indx[,1]
          tmp.index.control  <- indx[,2]

          tl  <- length(tmp.index.treated)
          index.treated <- vector(length=tl, mode="numeric")
          index.control <- vector(length=tl, mode="numeric")
          trt  <- Tr[tmp.index.treated]==1
          for (i in 1:tl)
            {
              if (trt[i]) {
                index.treated[i]  <- tmp.index.treated[i]
                index.control[i]  <- tmp.index.control[i]
              } else {
                index.treated[i]  <- tmp.index.control[i]
                index.control[i]  <- tmp.index.treated[i]
              }
            }
        } else if(estimand==2) {
          #"ATC"
          index.treated  <- indx[,2]
          index.control  <- indx[,1]
        }        

        mdata  <- list()
        mdata$Y  <- c(Y[index.treated],Y[index.control])
        mdata$Tr <- c(Tr[index.treated],Tr[index.control])
        mdata$X  <- rbind(X[index.treated,],X[index.control,])
        
        #naive standard errors
        mest  <- sum((Y[index.treated]-Y[index.control])*weights)/sum(weights)
        v1  <- Y[index.treated] - Y[index.control]
        varest  <- sum( ((v1-mest)^2)*weights)/(sum(weights)*sum(weights))
        se.naive  <- sqrt(varest)
        
        z  <- list(est=ret.use$est, se=ret.use$se, est.noadj=mest, se.naive=se.naive,
                   se.cond=ret.use$se.cond, 
                   mdata=mdata, em=ret.use$em,
                   index.treated=index.treated, index.control=index.control,
                   weights=weights, orig.nobs=orig.nobs, orig.wnobs=orig.wnobs,
                   orig.treated.nobs=orig.treated.nobs,
                   nobs=nobs, wnobs=sum(weights),
                   caliper=caliper, index.caliper.drop=index.drop, 
                   index.treated.nocaliper = index.treated.nocaliper,
                   index.control.nocaliper = index.control.nocaliper,
                   index.treated.indata=index.treated.indata,
                   index.control.indata=index.control.indata,
                   art.data=ret.use$art.data, aug.data=ret.use$aug.data)
      } else  {
        z  <- NA
      }

    class(z)  <- "Match"    
    return(z)
  } #end of Match

summary.Match  <- function(object, ..., full=FALSE, digits=5)
  {
    if(!is.list(object)) {
      warning("'Match' object contains no valid matches")
      return()
    }
    
    if (class(object) != "Match") {
      warning("Object not of class 'Match'")
      return()
    }

    cat("\n")
    cat("Estimate... ",format(object$est,digits=digits),"\n")
    cat("SE......... ",format(object$se,digits=digits),"\n")
    cat("T-stat..... ",format(object$est/object$se,digits=digits),"\n")
    cat("p.val...... ",format.pval((1-pnorm(abs(object$est/object$se)))*2,digits=digits),"\n")
    cat("\n")

    if(full)
      {
        cat("Est noAdj.. ",format(object$est.noadj,digits=digits),"\n")        
        cat("Naive SE... ",format(object$se.naive,digits=digits),"\n")
        cat("Naive T-st. ",format(object$est.noadj/object$se.naive,digits=digits),"\n")
        cat("Naive p.val ",format.pval((1-pnorm(abs(object$est.noadj/object$se.naive)))*2,digits=digits),"\n")
        cat("\n")
      }

    
    cat("Original number of observations (weighted)... ", object$orig.wnobs,"\n")
    cat("Original number of observations (unweighted). ", object$orig.nobs,"\n")
    cat("Original number of treated obs  (unweighted). ", object$orig.treated.nobs,"\n")    
    cat("Matched number of observations  (weighted)... ", signif(object$wnobs),"\n")
    cat("Matched number of observations  (unweighted). ", length(object$index.treated),"\n")

    cat("\n")
    if(object$caliper!=FALSE)
      {
        cat("Caliper (SDs)..............  ",object$caliper,"\n")
        cat("Matches dropped by caliper.. ",length(object$index.caliper.drop),"\n")
        cat("\n\n")
      } else {
        cat("\n")
      }
  } #end of summary.Match


Rmatch <- function(Y, Tr, X, Z, V, All, M, BiasAdj, Weight, Weight.matrix, Var.calc, weight,
                   SAMPLE, exact, ccc, cdd)
  {

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

# two slippage parameters that are used to determine whether distances are equal
# distances less than ccc or cdd are interpeted as zero.
# these are passed in.  ccc, cdd


# I. set up
# I.a. vector for which observations we want the average effect
# iot_t is the vector with weights in the average treatment effects
# iot_c is the vector of indicators for potential controls

    if (All==1)
      {
        iot.t <- weight;
        iot.c <- as.matrix(rep(1,length(Tr)))
      } else {
        iot.t <- Tr*weight;
        iot.c <- 1-Tr    
      }

# I.b. determine sample and covariate vector sizes
    N  <- nrow(X)
    Kx <- ncol(X)
    Kz <- ncol(Z)

# K covariates
# N observations
    Nt <- sum(Tr)
    Nc <- sum(1-Tr)
    on <- as.matrix(rep(1,N))

# I.c. normalize regressors to have mean zero and unit variance.
# If the standard deviation of a variable is zero, its normalization
# leads to a variable with all zeros.
# The matrix AA enables one to transform the user supplied weight matrix 
# to take account of this transformation.
# Mu_X and Sig_X keep track of the original mean and variances
    AA    <- diag(Kx)
    Mu.X  <- matrix(0, Kx, 1)
    Sig.X <- matrix(0, Kx, 1)

    for (k in 1:Kx)
      {
        Mu.X[k,1] <- sum(X[,k]*weight)/sum(weight)
        eps <- X[,k]-Mu.X[k,1]
        Sig.X[k,1] <- sqrt(sum(X[,k]*X[,k]*weight)/sum(weight)-Mu.X[k,1]^2)
        Sig.X[k,1] <- Sig.X[k,1]*sqrt(N/(N-1))
        X[,k]=eps/Sig.X[k,1]
        AA[k,k]=Sig.X[k,1]
      } #end of k loop

    Nv <- nrow(V)
    Mv <- ncol(V)
    Mu.V  <- matrix(0, Mv, 1)
    Sig.V <- matrix(0, Mv, 1)

    for (j in 1:Mv)
      {
        Mu.V[j,1]= ( t(V[,j])%*%weight ) /sum(weight)
        dv <- V[,j]-Mu.V[j,1]
        sv <- sum(V[,j]*V[,j]*weight)/sum(weight)-Mu.V[j,1]^2
        if (sv > 0)
          {
            sv <- sqrt(sv)
          } else {
            sv <- 0
          }
        sv <- sv * sqrt(N/(N-1))
        Sig.V[j,1] <- sv
      } #end of j loop

# I.d. define weight matrix for metric, taking into account normalization of
# regressors.
# If the eigenvalues of the regressors are too close to zero the Mahalanobis metric
# is not used and we revert back to the default inverse of variances.
    if (Weight==1)
      {
        Weight.matrix=diag(Kx)
      } else if (Weight==2) {
        if (min (eigen( t(X)%*%X/N, only.values=TRUE)$values) > 0.0000001)
          {
            Weight.matrix= solve(t(X)%*%X/N) 
          } else {
            Weight.matrix <- diag(Kx)
          }
      } else if (Weight==3)
        {
          Weight.matrix <- AA %*% Weight.matrix %*% AA
        }

    if (exact==1)
      {
        Weight.matrix <- cbind(Weight.matrix, matrix(0,nrow=Kx,ncol=Mv))
        Weight.matrix <- rbind(Weight.matrix, cbind(matrix(0,nrow=Mv,ncol=Kx),
                               1000*solve(diag(as.vector(Sig.V*Sig.V), nrow=length(Sig.V)))))
        Weight.matrix <- as.matrix(Weight.matrix)
        X <- cbind(X,V)
        Mu.X  <- rbind(Mu.X, matrix(0, nrow(Mu.V), 1))
        Sig.X <- rbind(Sig.X, matrix(1, nrow(Sig.V), 1))
      } #end of exact

    Nx <- nrow(X)
    Kx <- ncol(X)

    if ( min(eigen(Weight.matrix, only.values=TRUE)$values) < ccc )
      Weight.matrix <- Weight.matrix + diag(Kx)*ccc

    # I.fg. initialize matrices before looping through sample
    YCAUS  <- as.matrix(rep(0, N))
    SCAUS  <- as.matrix(rep(0, N))
    XCAUS  <- matrix(0, nrow=N, ncol=Kx)
    ZCAUS  <- matrix(0, nrow=N, ncol=Kz)
    Kcount <- as.matrix(rep(0, N))
    KKcount <- as.matrix(rep(0, N))
    MMi     <- as.matrix(rep(0, N))

# II. Loop through all observations that need to be matched.
    INN <- as.matrix(1:N)
    ww <- chol(Weight.matrix) # so that ww*ww=w.m
    TT <- as.matrix(1:N)

    # initialize some data objects
    DD <- NULL
    I  <- NULL
    IM <- NULL
    IT <- NULL
    IX <- NULL
    IZ <- NULL
    IY <- NULL
    W  <- NULL
    ADist <- NULL
    WWi <- NULL
    Xt <- NULL
    Zt <- NULL
    Yt <- NULL
    Xc <- NULL
    Zc <- NULL
    Yc <- NULL

    for (i in 1:N)
      {
        #treatment indicator for observation to be matched        
        TREATi <- Tr[i]   

        # proceed with all observations if All==1
        # but only with treated observations if All=0        
        if ( (TREATi==1 & All!=1) | All==1 )
          {
            # covariate value for observation to be matched                        
            xx <- t(as.matrix(X[i,]))

            # covariate value for observation to be matched            
            zz <- t(as.matrix(Z[i,]))

            # outcome value for observation to be matched                        
            yy <- Y[i]

            #JSS: check *
            foo <- as.matrix(rep(1, Nx))
            DX <- (X - foo %*% xx) %*% t(ww)

            if (Kx>1)
              {
                #JSS
                foo <- t(DX*DX)
                Dist <- as.matrix(apply(foo, 2, sum))
              } else {
                Dist <- as.matrix(DX*DX)
              } #end of Kx

            # Dist distance to observation to be matched
            # is N by 1 vector

            # set of potential matches (all observations with other treatment)
            # JSS, note:logical vector
            POTMAT <- Tr == (1-TREATi) 

            # X's for potential matches
            XPOT <- X[POTMAT,]
            DistPot <- Dist[POTMAT,1]
            TTPotMat <- TT[POTMAT,1]
            weightPot <- as.matrix(weight[POTMAT,1])

            # sorted distance of potential matches            
            S <- sort(DistPot)
            L <- order(DistPot)
            weightPot.sort <- weightPot[L,1]
            weightPot.sum <- cumsum(weightPot.sort)
            tt <- 1:(length(weightPot.sum))
            MMM <- min(tt[weightPot.sum>=M])

            # distance at last match
            Distmax <- S[MMM]

            # selection of actual matches 
            #logical index
            ACTMAT <- POTMAT & ( Dist <= (Distmax+ccc) )

            # distance to actual matches            
            ACTDIST <- as.matrix(Dist[ACTMAT,1])

            # counts how many times each observation is matched.
            Kcount <- Kcount + weight[i] * weight*ACTMAT/sum(ACTMAT*weight)
            
            KKcount <- KKcount+weight[i,1] * weight*weight*ACTMAT /
              (sum(ACTMAT*weight)*sum(ACTMAT*weight))

            # counts how many times each observation is matched.
            # counts number of matches for observation i
            # Unless there are ties this should equal M

            Mi <- sum(weight*ACTMAT)            
            MMi[i,1] <- Mi
            Wi <- as.matrix(weight[ACTMAT,1])

            # mean of Y's for actual matches
            ymat <- t(Y[ACTMAT,1]) %*% Wi/Mi

            # mean of X's for actual matches
            # mean of Z's for actual matches
            if (length(Wi)>1)
              {
                xmat <- t(t(X[ACTMAT,]) %*% Wi/Mi)
                zmat <- t(t(Z[ACTMAT,]) %*% Wi/Mi)
              } else {
                xmat <- t(X[ACTMAT,]) * as.real(Wi)/Mi
                zmat <- t(Z[ACTMAT,]) * as.real(Wi)/Mi
              }

            # estimate causal effect on y for observation i
            YCAUS[i,1] <- TREATi %*% (yy-ymat)+(1-TREATi) %*% (ymat-yy)

            # difference between x and actual matches for observation i
            XCAUS[i,] <- TREATi %*% (xx-xmat)+(1-TREATi) %*% (xmat-xx)

            ZCAUS[i,] <- TREATi %*% (zz-zmat)+(1-TREATi) %*% (zmat-zz)

            # collect results
            I  <- rbind(I, i * matrix(1, nrow=sum(ACTMAT), ncol=1))
            DD <- rbind(DD, ACTDIST)
            IM <- rbind(IM, as.matrix(INN[ACTMAT,1]))
            IT <- rbind(IT, TREATi * as.matrix(rep(1, sum(ACTMAT))))
            IX <- rbind(IX, as.matrix(rep(1, sum(ACTMAT))) %*% xx)
            IZ <- rbind(IZ, as.matrix(rep(1, sum(ACTMAT))) %*% zz)
            IY <- rbind(IY, as.matrix(rep(1, sum(ACTMAT))) * yy)

            # weight for matches
            W <- as.matrix(c(W, weight[i,1] * Wi/Mi))
            ADist <- as.matrix(c(ADist, ACTDIST))
            WWi   <- as.matrix(c(WWi, Wi))

            if (TREATi==1)
              {

                if (ncol(X) > 1)
                  {
                    # covariates for treated
                    Xt <- rbind(Xt, as.matrix(rep(1, sum(ACTMAT))) %*% xx)                    
                    # covariate for matches
                    Xc <- rbind(Xc, X[ACTMAT,])
                    
                  } else {
                    # covariates for treated
                    Xt <- as.matrix(c(Xt, as.matrix(rep(1, sum(ACTMAT))) %*% xx))
                    # covariate for matches
                    Xc <- as.matrix(c(Xc, X[ACTMAT,]))
                  }

                if (ncol(Z) > 1)
                  {
                    # covariates for treated                
                    Zt <- rbind(Zt, as.matrix(rep(1, sum(ACTMAT))) %*% zz)
                    # covariate for matches
                    Zc <- rbind(Zc, Z[ACTMAT,])
                  } else {
                    # covariates for treated                
                    Zt <- as.matrix(c(Zt, as.matrix(rep(1, sum(ACTMAT))) %*% zz))
                    # covariate for matches                
                    Zc <- as.matrix(c(Zc, Z[ACTMAT,]))
                  }

                # outcome for treated                              
                Yt <- as.matrix(c(Yt, yy * as.matrix(rep(1, sum(ACTMAT))) ))
                # outcome for matches
                Yc <- as.matrix(c(Yc, Y[ACTMAT,1]))
              } else {

                if (ncol(X) > 1)
                  {
                    # covariates for controls
                    Xc <- rbind(Xc, as.matrix(rep(1, sum(ACTMAT))) %*% xx)
                    # covariate for matches
                    Xt <- rbind(Xt, X[ACTMAT,])
                  } else {
                    # covariates for controls                
                    Xc <- as.matrix(c(Xc, as.matrix(rep(1, sum(ACTMAT))) %*% xx))
                    # covariate for matches
                    Xt <- as.matrix(c(Xt, X[ACTMAT,]))
                  }

                if (ncol(Z) > 1)
                  {
                    # covariates for controls                
                    Zc <- rbind(Zc, as.matrix(rep(1, sum(ACTMAT))) %*% zz)
                    # covariate for matches                
                    Zt <- rbind(Zt, Z[ACTMAT,])
                  } else {
                    # covariates for controls                
                    Zc <- as.matrix(c(Zc, as.matrix(rep(1, sum(ACTMAT))) %*% zz))
                    # covariate for matches                
                    Zt <- as.matrix(c(Zt, Z[ACTMAT,]))
                  }
                # outcome for controls                               
                Yc <- as.matrix(c(Yc, as.matrix(rep(1, sum(ACTMAT))) * yy))
                # outcome for matches
                Yt <- as.matrix(c(Yt, Y[ACTMAT,1]))

              } #end of TREATi
          } #end of if
      } #i loop

    # transform matched covariates back for artificial data set
    Xt.u <- Xt
    Xc.u <- Xc
    IX.u <- IX
    for (k in 1:Kx)
      {
        Xt.u[,k] <- Mu.X[k,1]+Sig.X[k,1] * Xt.u[,k]
        Xc.u[,k] <- Mu.X[k,1]+Sig.X[k,1] * Xc.u[,k]
        IX.u[,k] <- Mu.X[k,1]+Sig.X[k,1] * IX.u[,k]
      }
    if (All!=1)
      {
        I  <- as.matrix(I[IT==1,])
        IM <- as.matrix(IM[IT==1,])
        IT <- as.matrix(IT[IT==1,])
        IY <- as.matrix(IY[IT==1,])
        Yc <- as.matrix(Yc[IT==1,])
        Yt <- as.matrix(Yt[IT==1,])
        W  <- as.matrix(W[IT==1,])
        ADist <- as.matrix(ADist[IT==1,])
        WWi   <- as.matrix(WWi[IT==1,])
        IX.u  <- as.matrix(IX.u[IT==1,])
        Xc.u  <- as.matrix(Xc.u[IT==1,])
        Xt.u  <- as.matrix(Xt.u[IT==1,])
        Xc    <- as.matrix(Xc[IT==1,])
        Xt <- as.matrix(Xt[IT==1,])
        Zc <- as.matrix(Zc[IT==1,])
        Zt <- as.matrix(Zt[IT==1,])
        IZ <- as.matrix(IZ[IT==1,])
      } #end of if

    # III. Regression of outcome on covariates for matches
    if (All==1)
      {
        # number of observations
        NNt <- nrow(Z)
        # add intercept        
        ZZt <- cbind(rep(1, NNt), Z)
        # number of covariates
        Nx <- nrow(ZZt)
        Kx <- ncol(ZZt)
        xw <- ZZt*(sqrt(Tr*Kcount) %*% t(as.matrix(rep(1,Kx))))

        foo <- min(eigen(t(xw)%*%xw, only.values=TRUE)$values)
        foo <- as.real(foo<=cdd)
        foo2 <- apply(xw, 2, sd)

        wout <- solve( t(xw) %*% xw + diag(Kx) * cdd * (foo) * max(foo2)) %*%
          (t(xw) %*% (Y*sqrt(Tr*Kcount)))

        NW <- nrow(wout)
        KW <- ncol(wout)
        Alphat <- wout[2:NW,1]
      } else {
        Alphat <- matrix(0, nrow=Kz, ncol=1)
      } #end if ALL

    # III.b.  Controls
    NNc <- nrow(Z)
    ZZc <- cbind(matrix(1, nrow=NNc, ncol=1),Z)
    Nx <- nrow(ZZc)
    Kx <- ncol(ZZc)

    xw <- ZZc*(sqrt((1-Tr)*Kcount) %*% matrix(1, nrow=1, ncol=Kx))

    foo <- min(eigen(t(xw)%*%xw, only.values=TRUE)$values)
    foo <- as.real(foo<=cdd)
    foo2 <- apply(xw, 2, sd)

    wout <- solve( t(xw) %*% xw + diag(Kx) * cdd * (foo) * max(foo2)) %*%
      (t(xw) %*% (Y*sqrt((1-Tr)*Kcount)))
    
    NW <- nrow(wout)
    KW <- ncol(wout)
    Alphac <- as.matrix(wout[2:NW,1])

    Alpha <- cbind(Alphat,Alphac)

    # III.c. adjust matched outcomes using regression adjustment for bias adjusted matching estimator

    SCAUS <- YCAUS-Tr*(ZCAUS %*% Alphac)-(1-Tr)*(ZCAUS %*% Alphat)
    # adjusted treated outcome
    Yc.adj <- Yc+BiasAdj * (IZ-Zc) %*% Alphac
    # adjusted control outcome
    Yt.adj <- Yt+BiasAdj*(IZ-Zt) %*% Alphat
    Yt.adj <- Yt+BiasAdj*(IZ-Zt) %*% Alphat
    Tau.i <- Yt.adj - Yc.adj

    art.data <- cbind(I,IM,IT,DD,IY,Yc,Yt,W,WWi,ADist,IX.u,Xc.u,Xt.u,
                      Yc.adj,Yt.adj,Tau.i)

    # III. If conditional variance is needed, initialize variance vector
    # and loop through all observations

    Nx <- nrow(X)
    Kx <- ncol(X)

    ww <- chol(Weight.matrix)
    NN <- as.matrix(1:N)
    if (Var.calc>0)
      {
        Sig <- matrix(0, nrow=N, ncol=1)
        # overall standard deviation of outcome
        # std <- sd(Y)
        for (i in 1:N)
          {
            # treatment indicator observation to be matched
            TREATi <- Tr[i,1]
            # covariate value for observation to be matched
            xx <- X[i,]
            # potential matches are all observations with the same treatment value
            POTMAT <- (Tr==TREATi)
            POTMAT[i,1] <- 0
            weightPOT <- as.matrix(weight[POTMAT==1,1])
            DX <- (X - matrix(1, Nx,1) %*% xx) %*% t(ww)
            if (Kx>1)
              {
                foo <- apply(t(DX*DX), 2, sum)
                Dist <- as.matrix(foo)
              } else {
                Dist <- DX*DX
              }

            # distance to observation to be matched

            # Distance vector only for potential matches
            DistPot <- Dist[POTMAT==1,1]
            # sorted distance of potential matches
            S <- sort(DistPot)
            L <- order(DistPot)
            weightPOT.sort <- weightPOT[L,1]
            weightPOT.sum <- cumsum(weightPOT.sort)
            tt <-  1:(length(weightPOT.sum))
            MMM <- min(tt[weightPOT.sum >= Var.calc])
            MMM <- min(MMM,length(S))
            Distmax=S[MMM]

            # distance of Var_calc-th closest match
            ACTMAT <- (POTMAT==1) & (Dist<= (Distmax+ccc))

            # indicator for actual matches, that is all potential
            # matches closer than, or as close as the Var_calc-th
            # closest

            Yactmat <- as.matrix(c(Y[i,1], Y[ACTMAT,1]))
            weightactmat <- as.matrix(c(weight[i,1], weight[ACTMAT,1]))
            fm <- t(Yactmat) %*% weightactmat/sum(weightactmat)
            sm <- sum(Yactmat*Yactmat*weightactmat)/sum(weightactmat)
            sigsig <- (sm-fm %*% fm)*sum(weightactmat)/(sum(weightactmat)-1)
            
            # standard deviation of actual matches
            Sig[i,1] <- sqrt(sigsig)
          }# end of i loop
        #variance estimate 
        Sigs <- Sig*Sig
      } #end of var.calc > 0

    # matching estimator
    est <- t(W) %*% Tau.i/sum(W)
    est.t <- sum((iot.t*Tr+iot.c*Kcount*Tr)*Y)/sum(iot.t*Tr+iot.c*Kcount*Tr)
    est.c <- sum((iot.t*(1-Tr)+iot.c*Kcount*(1-Tr))*Y)/sum(iot.t*(1-Tr)+iot.c*Kcount*(1-Tr))

    if (Var.calc==0)
      {
        eps <- Tau.i - as.real(est)
        eps.sq <- eps*eps
        Sigs <- 0.5 * matrix(1, N, 1) %*% (t(eps.sq) %*% W)/sum(W)
        sss <- sqrt(Sigs[1,1])
      } #end of Var.calc==0

    SN <- sum(iot.t)
    var.sample <- sum((Sigs*(iot.t+iot.c*Kcount)*(iot.t+iot.c*Kcount))/(SN*SN))

    if (All==1)
      {
        var.pop <- sum((Sigs*(iot.c*Kcount*Kcount+2*iot.c*Kcount-iot.c*KKcount))/(SN*SN))
      } else {
        var.pop=sum((Sigs*(iot.c*Kcount*Kcount-iot.c*KKcount))/(SN*SN))
      }

    if (BiasAdj==1)
      {
        dvar.pop <- sum(iot.t*(SCAUS-as.real(est))*(SCAUS-as.real(est)))/(SN*SN)
      } else {
        dvar.pop <- sum(iot.t*(YCAUS-as.real(est))*(YCAUS-as.real(est)))/(SN*SN)
      }

    var.pop <- var.pop + dvar.pop

    if (SAMPLE==1)
      {
        var <- var.sample
      } else {
        var <- max(var.sample, var.pop)
        var <- var.pop
      }

    var.cond <- max(var.sample,var.pop)-var.sample
    se <- sqrt(var)
    se.cond <- sqrt(var.cond)

    Sig <- sqrt(Sigs)
    aug.data <- cbind(Y,Tr,X,Z,Kcount,Sig,weight)

    if (All==2)
      est <- -1*est

    em <- as.matrix(c(0,0))

    if (exact==1)
      {
        Vt.u <- Xt.u[,(Kx-Mv+1):Kx]
        Vc.u <- Xc.u[,(Kx-Mv+1):Kx]
        Vdif <- abs(Vt.u-Vc.u)

        if (Mv>1)
          Vdif <- as.matrix(apply(t(Vdif), 2, sum))

        em[1,1] <- length(Vdif)
        em[2,1] <- sum(Vdif>0.000001)
      }#end of exact==1

    return(list(est=est, se=se, se.cond=se.cond, em=em, W=W,
                art.data=art.data, aug.data=aug.data))
  }# end of Rmath


#
# See Rosenbaum and Rubin (1985) and Smith and Todd Rejoinder (JofEconometrics) p.9
#
sdiff  <- function(Tr, Co, weights=rep(1,length(Co)))
  {
    equal  <- FALSE
    if (length(Tr)==length(Co))
      equal  <- TRUE

    if (!equal)
      {
        diff  <- mean(Tr)-mean(Co)
        sdiff <- 100*diff/sqrt( (var(Tr)+var(Co))/2 )
      } else{
        diff  <- sum( (Tr-Co)*weights ) /sum(weights)
        mean.Tr <- sum(Tr*weights)/sum(weights)
        mean.Co <- sum(Co*weights)/sum(weights)
        var.Tr  <- sum( ( (Tr - mean.Tr)^2 )*weights)/sum(weights)
        var.Co  <- sum( ( (Co - mean.Co)^2 )*weights)/sum(weights)
        sdiff  <- 100*diff/sqrt( (var.Tr+var.Co)/2 )
      }

    return(sdiff)
  }

# function runs sdiff and t.test
#
balanceUV  <- function(Tr, Co, weights=rep(1,length(Co)), exact=FALSE, ks=FALSE, nboots=1000)
  {
    try.McNemar <- FALSE
    equal  <- FALSE
    if (length(Tr)==length(Co))
      equal  <- TRUE

    ks.out  <- NULL

    categorical  <- is.factor(Tr)

    if (categorical & equal)
      {
        x <- factor(Tr)
        y <- factor(Co)
        r <- nlevels(x)
        if ((r < 2) || (nlevels(y) != r))
          {
            cat("reverting to wilcox.test: x and y must have the same number of levels (minimum 2)\n")
            try.McNemar  <- FALSE
          } else {
            efoo  <- sum(levels(x)==levels(y))
            if (efoo != r)
              {
                cat("reverting to wilcox.test: Level sets of factors are different\n")
                try.McNemar  <- FALSE            
              }
          }
      }    

    if (equal & categorical & try.McNemar)
      {
        use.McNemar  <- TRUE
      } else {
        use.McNemar  <- FALSE
        Tr  <- as.real(Tr)
        Co  <- as.real(Co)
      }

    if (!equal)
      {
        sbalance  <- sdiff(Tr, Co, weights)
        
        mean.Tr  <- mean(Tr);
        mean.Co  <- mean(Co);
        var.Tr  <- var(Tr)
        var.Co  <- var(Co)
        var.ratio  <- var.Tr/var.Co

        wc  <- wilcox.test(Tr,Co, exact=exact)
        tt  <- t.test(Tr, Co)

        if (ks)
          ks.out <- ks.boot(Tr, Co,nboots=nboots)

        pdiscordant  <- NA
      } else {

        if (categorical & use.McNemar)
          {
            sbalance  <- NULL
            
            x  <- as.real(Tr)
            y  <- as.real(Co)
            mean.Tr  <- sum(x*weights)/sum(weights);
            mean.Co  <- sum(y*weights)/sum(weights);
            var.Tr  <- sum( ( (x - mean.Tr)^2 )*weights)/sum(weights);
            var.Co  <- sum( ( (y - mean.Co)^2 )*weights)/sum(weights);
            var.ratio  <- var.Tr/var.Co


            mc  <- McNemar2(Tr, Co, weights=weights)
            wc  <- wilcox.test( (as.real(Tr)-as.real(Co))*weights, exact=exact )

            if (ks)
              ks.out  <- ks.boot(Tr, Co,nboots=nboots)            

          } else {
            sbalance  <- sdiff(Tr, Co, weights)
            
            mean.Tr  <- sum(Tr*weights)/sum(weights);
            mean.Co  <- sum(Co*weights)/sum(weights);
            var.Tr  <- sum( ( (Tr - mean.Tr)^2 )*weights)/sum(weights);
            var.Co  <- sum( ( (Co - mean.Co)^2 )*weights)/sum(weights);
            var.ratio  <- var.Tr/var.Co

            wc  <- wilcox.test( (Tr-Co)*weights, exact=exact )

#            tt  <- t.test((Tr-Co)*weights)
            tt  <- Mt.test(Tr, Co, weights)

            if (ks)
              ks.out  <- ks.boot(Tr, Co, nboots=nboots)
          }        
      }


    if (use.McNemar)
      {
        ret  <- list(sdiff=sbalance,mean.Tr=mean.Tr,mean.Co=mean.Co,
                     var.Tr=var.Tr,var.Co=var.Co, p.value=mc$p.value,
                     var.ratio=var.ratio,
                     pdiscordant=mc$pdiscordant, ks=ks.out,
                     wc=NULL, mc=mc, tt=NULL)        
      } else {
        ret  <- list(sdiff=sbalance,mean.Tr=mean.Tr,mean.Co=mean.Co,
                     var.Tr=var.Tr,var.Co=var.Co, p.value=tt$p.value,
                     var.ratio=var.ratio,
                     pdiscordant=NULL, ks=ks.out,
                     wc=wc, mc=NULL, tt=tt)
      }

    class(ret) <-  "balanceUV"
    return(ret)
  } #balanceUV

summary.balanceUV  <- function(object, ..., digits=5)
  {

    if (class(object) != "balanceUV") {
      warning("Object not of class 'balanceUV'")
      return(NULL)
    }

    use.McNemar <- FALSE
    if (is.list(object$mc))
      use.McNemar <- TRUE

    if (use.McNemar)
      {
        cat("mean treatment........", format(object$mean.Tr, digits=digits),"\n")
        cat("mean control..........", format(object$mean.Co, digits=digits),"\n")            
#       cat("Wilcoxon pval.........", format.pval(wc$p.value,digits=digits), "\n")
        cat("var ratio (Tr/Co).....", format(object$var.ratio, digits=digits),"\n")
        cat("Discordant prop.......", format(object$pdiscordant,digits=digits),"\n")
        cat("McNemar pval..........", format.pval(object$p.value,digits=digits), "\n")
        if (!is.null(object$ks$ks.boot.pvalue))
          {
            cat("KS Bootstrap p-value..", format.pval(object$ks$ks.boot.pvalue, digits=digits), "\n")
            cat("KS Statistic..........", format(object$ks$ks$statistic, digits=digits), "\n")
            cat("KS Naive p-value......", format(object$ks$ks$p.value, digits=digits), "\n")            
          }
      } else {
        cat("mean treatment........", format(object$mean.Tr, digits=digits),"\n")
        cat("mean control..........", format(object$mean.Co, digits=digits),"\n")
#       cat("sdiff.................", format(object$sdiff, digits=digits),"\n")            
#       cat("Wilcoxon pval.........", format.pval(wc$p.value,digits=digits), "\n")
        cat("var ratio (Tr/Co).....", format(object$var.ratio, digits=digits),"\n")
        cat("T-test pval...........", format.pval(object$tt$p.value,digits=digits), "\n")            
        if (!is.null(object$ks$ks.boot.pvalue))
          {
            cat("KS Bootstrap p-value..", format.pval(object$ks$ks.boot.pvalue, digits=digits), "\n")
            cat("KS Statistic..........", format(object$ks$ks$statistic, digits=digits), "\n")
            cat("KS Naive p-value......", format(object$ks$ks$p.value, digits=digits), "\n")            
          }        
      }
    
    cat("\n")        
  } #end of summary.balanceUV




McNemar  <- function(Tr, Co, weights=rep(1,length(Tr)))
{
#SEE MCNEMAR2 FOR GENERAL FUNCTION  
#McNemar's test
#mcnemar.test(x, y = NULL, correct = TRUE)
#McNemar's test of symmetry is appropriate for parallel measures from
#matched cases as well as for repeated measures on a single set of
#cases. In this case, McNemar's test of symmetry is equivalent to
#Cochran's Q.
#
#Siegel, S. Nonparametric Methods for the Behavioral Sciences. New
#York: McGraw-Hill, 1956. p. 63 (when both variables are two-point
#scales, McNemar's test of symmetry and McNemar's test for the
#significance of changes are equivalent); Bowker, A. H., A test for
#symmetry in contingency tables. Journal of the American Statistical
#Association 43 (1948): 572-574.

#page 148ff Argesti
  
  obs  <- sum(weights);
  p1.  <- sum(Tr*weights)/obs;
  p.1  <- sum(Co*weights)/obs;
  p1  <- p1.;
  p2  <- p.1;
  mc.table.obs  <-  as.matrix(table(Tr, Co));
  mc.table  <-  mc.table.obs/obs;
  p11  <- mc.table[1,1];
  p12  <- mc.table[1,2];
  p21  <- mc.table[2,1];
  p22  <- mc.table[2,2];
  
  #Vd  <- (p1.*(1-p1.) + p.1*(1-p.1) - 2*(p11*p22-p12*p21))/obs;
  #Sd  <- sqrt(Vd);
  #d   <- p1-p2;

  #
  #cat("Mcnemar 10.2 test: ", d/Sd, (1-pnorm(abs(d/Sd)))*2,"\n");

  n11.indx  <- Tr==TRUE  & Co==TRUE
  n12.indx  <- Tr==TRUE  & Co==FALSE
  n21.indx  <- Tr==FALSE & Co==TRUE
  n22.indx  <- Tr==FALSE & Co==FALSE
    
  n11  <- sum(n11.indx*weights)
  n12  <- sum(n12.indx*weights)
  n21  <- sum(n21.indx*weights)
  n22  <- sum(n22.indx*weights)
  
  mc  <- (n12-n21)/sqrt(n12 + n21);
  pdiscordant= (n12+n21)/obs

  #this matches the results of mcnemar.test without the continuity correction
  #cat("Mcnemar 10.3 test: ", mc, ,"\n");
  
  # mcnemar  <- mcnemar.test(Tr, Co,correct = FALSE);
  p.value  <- (1-pnorm(abs(mc)))*2;
  return(list(statistic=mc,p.value=p.value,pdiscordant=pdiscordant))
}


McNemar2 <- function (Tr, Co, correct = TRUE, weights=rep(1,length(Tr)))
{
  x  <- Tr
  y  <- Co
  if (is.matrix(x)) {
    stop("this version of McNemar cannot handle x being a matrix")
  } else {
    if (is.null(y)) 
      stop("if x is not a matrix, y must be given")
    if (length(x) != length(y)) 
      stop("x and y must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- factor(x[OK])
    y <- factor(y[OK])
    r <- nlevels(x)
    if ((r < 2) || (nlevels(y) != r))
      {
        stop("x and y must have the same number of levels (minimum 2)")
      }
  }
  tx <- table(x, y)
  facs  <- levels(x)
  txw  <- tx
  for(i in 1:r)
    {
      for(j in 1:r)
        {
          indx  <- x==facs[i] & y==facs[j]
          txw[i,j]  <- sum(weights[indx]);
        }
    }
  pdiscordant  <- sum( ( (x!=y)*weights )/sum(weights) )
  x  <- txw
  
  PARAMETER <- r * (r - 1)/2
  METHOD <- "McNemar's Chi-squared test"
  if (correct && (r == 2) && any(x - t(x))) {
    y <- (abs(x - t(x)) - 1)
    METHOD <- paste(METHOD, "with continuity correction")
  } else y <- x - t(x)
  x <- x + t(x)
  STATISTIC <- sum(y[upper.tri(x)]^2/x[upper.tri(x)])
  PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
  names(STATISTIC) <- "McNemar's chi-squared"
  names(PARAMETER) <- "df"

  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
               p.value = PVAL, method = METHOD, data.name = DNAME,
               pdiscordant=pdiscordant)
  class(RVAL) <- "htest"
  return(RVAL)
}

ks<-function(x,y,w=F,sig=T){
#  Compute the Kolmogorov-Smirnov test statistic
#
# Code for the Kolmogorov-Smirnov test is adopted from the Splus code
# written by Rand R. Wilcox for his book titled "Introduction to
# Robust Estimation and Hypothesis Testing." Academic Press, 1997.
#
#  Also see 
#@book( knuth1998,
#  author=       {Knuth, Donald E.},
#  title=        {The Art of Computer Programming, Vol. 2: Seminumerical Algorithms},
#  edition=      "3rd",
#  publisher=    "Addison-Wesley", address= "Reading: MA", 
#  year=         1998
#)
# and Wilcox 1997 
#
#  w=T computes the weighted version instead. 
#  sig=T indicates that the exact significance level is to be computed.
#  If there are ties, the reported significance level is exact when
#  using the unweighted test, but for the weighted test the reported
#  level is too high.
#
#  This function uses the functions ecdf, kstiesig, kssig and kswsig
#
#  This function returns the value of the test statistic, the approximate .05
#  critical value, and the exact significance level if sig=T.
#
#  Missing values are automatically removed
#
x<-x[!is.na(x)]
y<-y[!is.na(y)]
w<-as.logical(w)
sig<-as.logical(sig)
tie<-logical(1)
tie<-F
siglevel<-NA
z<-sort(c(x,y))  # Pool and sort the observations
for (i in 2:length(z))if(z[i-1]==z[i])tie<-T #check for ties
v<-1   # Initializes v
for (i in 1:length(z))v[i]<-abs(ecdf(x,z[i])-ecdf(y,z[i]))
ks<-max(v)
crit<-1.36*sqrt((length(x)+length(y))/(length(x)*length(y))) # Approximate
#                                                       .05 critical value
if(!w && sig && !tie)siglevel<-kssig(length(x),length(y),ks)
if(!w && sig && tie)siglevel<-kstiesig(x,y,ks)
if(w){
crit<-(max(length(x),length(y))-5)*.48/95+2.58+abs(length(x)-length(y))*.44/95
if(length(x)>100 || length(y)>100){
print("When either sample size is greater than 100,")
print("the approximate critical value can be inaccurate.")
print("It is recommended that the exact significance level be computed.")
}
for (i in 1:length(z)){
temp<-(length(x)*ecdf(x,z[i])+length(y)*ecdf(y,z[i]))/length(z)
temp<-temp*(1.-temp)
v[i]<-v[i]/sqrt(temp)
}
v<-v[!is.na(v)]
ks<-max(v)*sqrt(length(x)*length(y)/length(z))
if(sig)siglevel<-kswsig(length(x),length(y),ks)
if(tie && sig){
print("Ties were detected. The reported significance level")
print("of the weighted Kolmogorov-Smirnov test statistic is not exact.")
}}
#round off siglevel in a nicer way
if(is.real(ks) & is.real(crit) & !is.na(ks) & !is.na(crit))
  {
    if (is.na(siglevel) & ks < crit)
      {
        siglevel  <- 0.99999837212332
      }
    if (is.real(siglevel) & !is.na(siglevel))
      {
        if (siglevel < 0)
          siglevel  <- 0
      }
  }

list(test=ks,critval=crit,pval=siglevel)
}


kssig<-function(m,n,val){
#
#    Compute significance level of the  Kolmogorov-Smirnov test statistic
#    m=sample size of first group
#    n=sample size of second group
#    val=observed value of test statistic
#
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 0:m){
    for (j in 0:n)cmat[i+1,j+1]<-abs(i/m-j/n)
  }
  cmat<-ifelse(cmat<=val,1e0,0e0)
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kssig<-1.-umat[m+1,n+1]/exp(term)
  return(kssig)
}

kstiesig<-function(x,y,val){
#
#    Compute significance level of the  Kolmogorov-Smirnov test statistic
#    for the data in x and y.
#    This function allows ties among the  values.
#    val=observed value of test statistic
#
m<-length(x)
n<-length(y)
z<-c(x,y)
z<-sort(z)
cmat<-matrix(0,m+1,n+1)
umat<-matrix(0,m+1,n+1)
for (i in 0:m){
for (j in 0:n){
if(abs(i/m-j/n)<=val)cmat[i+1,j+1]<-1e0
k<-i+j
if(k > 0 && k<length(z) && z[k]==z[k+1])cmat[i+1,j+1]<-1
}
}
for (i in 0:m){
for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
}
term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
kstiesig<-1.-umat[m+1,n+1]/exp(term)
kstiesig
}

ecdf<-function(x,val){
#  compute empirical cdf for data in x evaluated at val
#  That is, estimate P(X <= val)
#
ecdf<-length(x[x<=val])/length(x)
ecdf
}

kswsig<-function(m,n,val){
#
#    Compute significance level of the weighted
#    Kolmogorov-Smirnov test statistic
#
#    m=sample size of first group
#    n=sample size of second group
#    val=observed value of test statistic
#
mpn<-m+n
cmat<-matrix(0,m+1,n+1)
umat<-matrix(0,m+1,n+1)
for (i in 1:m-1){
for (j in 1:n-1)cmat[i+1,j+1]<-abs(i/m-j/n)*sqrt(m*n/((i+j)*(1-(i+j)/mpn)))
}
cmat<-ifelse(cmat<=val,1,0)
for (i in 0:m){
for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
}
term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
kswsig<-1.-umat[m+1,n+1]/exp(term)
kswsig
}

Mks.test <- function (x, y, ..., alternative = c("two.sided", "less", "greater"), 
                      exact = NULL, MC=FALSE) 
{
  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1) 
    stop("Not enough x data")
  PVAL <- NULL
  if (is.numeric(y)) {
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    y <- y[!is.na(y)]
    n.x <- as.double(n)
    n.y <- length(y)
    if (n.y < 1) 
      stop("Not enough y data")
    if (is.null(exact)) 
      exact <- (n.x * n.y < 10000)
    METHOD <- "Two-sample Kolmogorov-Smirnov test"
    TIES <- FALSE
    n <- n.x * n.y/(n.x + n.y)
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    if (length(unique(w)) < (n.x + n.y)) {
      if (!MC)
        warning("Use the Monte Carlo KS test because the regular test cannot compute correct p-values with ties")
      z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
      TIES <- TRUE
    }
    STATISTIC <- switch(alternative, two.sided = max(abs(z)), 
                        greater = max(z), less = -min(z))
    if (exact && alternative == "two.sided" && !TIES) 
      PVAL <- 1 - .C("psmirnov2x", p = as.double(STATISTIC), 
                     as.integer(n.x), as.integer(n.y), PACKAGE = "stats")$p
  }
  else {
    if (is.character(y)) 
      y <- get(y, mode = "function")
    if (mode(y) != "function") 
      stop("y must be numeric or a string naming a valid function")
    METHOD <- "One-sample Kolmogorov-Smirnov test"
    if (length(unique(x)) < n)
      {
        if (!MC)
          warning("Use the Monte Carlo KS test because the regular test cannot compute correct p-values with ties")
      }
    x <- y(sort(x), ...) - (0:(n - 1))/n
    STATISTIC <- switch(alternative, two.sided = max(c(x, 
                                       1/n - x)), greater = max(1/n - x), less = max(x))
  }
  names(STATISTIC) <- switch(alternative, two.sided = "D", 
                             greater = "D^+", less = "D^-")
  pkstwo <- function(x, tol = 1e-06) {
    if (is.numeric(x)) 
      x <- as.vector(x)
    else stop("Argument x must be numeric")
    p <- rep(0, length(x))
    p[is.na(x)] <- NA
    IND <- which(!is.na(x) & (x > 0))
    if (length(IND) > 0) {
      p[IND] <- .C("pkstwo", as.integer(length(x)), p = as.double(x[IND]), 
                   as.double(tol), PACKAGE = "stats")$p
    }
    return(p)
  }
  if (is.null(PVAL)) {
    PVAL <- ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) * 
                     STATISTIC), exp(-2 * n * STATISTIC^2))
  }
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = alternative, 
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
} #end of Mks.test

Mt.test  <- function(Tr, Co, weights)
  {
    v1  <- Tr-Co
    estimate  <- sum(v1*weights)/sum(weights)
    var1  <- sum( ((v1-estimate)^2)*weights )/( sum(weights)*sum(weights) )

    statistic  <- estimate/sqrt(var1)
    parameter  <- Inf
#    p.value    <- (1-pnorm(abs(statistic)))*2
    p.value    <- (1-pt(abs(statistic), df=sum(weights)-1))*2

    z  <- list(statistic=statistic, parameter=parameter, p.value=p.value,
               estimate=estimate)
    return(z)
  } #end of Mt.test



MatchBalance <- function(formul, data=NULL, match.out=NULL, ks=FALSE, mv=FALSE, 
                         nboots=1000, nmc=nboots, 
                         maxit=1000, weights=rep(1,nrow(data)),
                         digits=5, verbose=1, ...)
  {

    if(!is.list(match.out) & !is.null(match.out)) {
      warning("'Match' object contains no valid matches")
      match.out  <- NULL
    }

    if ( (class(match.out) != "Match") & (!is.null(match.out)) ) {
      warning("Object not of class 'Match'")
      match.out  <- NULL
    }

    if (is.null(data))
      {
        datain <- NULL
        xdata <- as.data.frame(get.xdata(formul,datafr=environment(formul)))
        ydata <- as.data.frame(get.ydata(formul,datafr=environment(formul)))
        Tr    <- ydata
        data  <- cbind(ydata,xdata)        
      } else {
        data  <- as.data.frame(data)
        datain <- data
        xdata  <- as.data.frame(get.xdata(formul, data))
        Tr  <- get.ydata(formul, data)        
      }

    if (sum(is.na(data)!=0))
      stop("MatchBalance: NAs found in data input")

    nvars  <- ncol(xdata)
    names.xdata  <- names(xdata)

    
    if (verbose > 0)
      {
        
        findx  <- 1
        if (sum(xdata[,1]==rep(1,nrow(xdata)))==nrow(xdata))
          {
            findx  <- 2
          }

        for( i in findx:nvars)
          {
            cat("\n***** (V",i-findx+1,") ", names.xdata[i]," *****\n",sep="")

            ks.do  <- FALSE
            is.dummy  <- length(unique( xdata[,i] )) < 3
            if (ks & !is.dummy)
              ks.do  <- TRUE

            cat("before matching:\n")
            foo  <-  balanceUV(xdata[,i][Tr==1],xdata[,i][Tr==0], ks=ks.do, nboots=nboots)
            summary(foo, digits=digits)
            
            if (!is.null(match.out))
              {
                cat("after  matching:\n")
                foo  <- balanceUV(xdata[,i][match.out$index.treated],
                                  xdata[,i][match.out$index.control],
                                  weights=match.out$weights, ks=ks.do, nboots=nboots)
                summary(foo, digits=digits)                
              }
          }
        
        if (mv)  {
          cat("...estimating Kolmogorov-Smirnov tests...\n")
          if (nboots > 0)
            cat("  ",nboots,"bootstraps are being run (see the 'nboots' option) \n")
          if (nmc > 0)
            cat("  ",nmc,"Monte Carlos are being run for the KS test in each bootstrap (see the 'nmc' option)\n")
          if ( (nboots < 500) | (nmc < 500) )
            {
              cat("\n  ","For publication quality p-values it is recommended that 'nboots' and 'nmc'\n")
              cat("  ","be set equal to at least 500 (preferably 1000).\n")              
            }
          cat("\n")
          
          ml  <- balanceMV(formul=formul, data=datain, match.out=match.out,
                           maxit=maxit, weights=weights, nboots=nboots, nmc=nmc,
                           verbose=verbose, ...)
          summary.balanceMV(ml, digits=digits)
        }  else  {
          ml  <- NULL
        }
      } else {
        for( i in 2:nvars)
          {
            ks.do  <- FALSE
            is.dummy  <- length(unique( xdata[,i] )) < 3
            if (ks & !is.dummy)
              ks.do  <- TRUE

            foo  <-  balanceUV(xdata[,i][Tr==1],xdata[,i][Tr==0],
                               weights=weights, ks=ks.do, nboots=nboots)
            
            if (!is.null(match.out))
              {
                foo  <- balanceUV(xdata[,i][match.out$index.treated],
                                  xdata[,i][match.out$index.control],
                                  weights=match.out$weights, ks=ks.do, nboots=nboots)
              }
          }
        
        if (mv)  {
          ml  <- balanceMV(formul=formul, data=datain, match.out=match.out,
                           maxit=maxit, weights=weights, nboots=nboots, nmc=nmc,
                           verbose=verbose, ...)
        }  else  {
          ml  <- NULL
        }        
      }#verbose end
    
    return(invisible(list(mv=ml,uv=foo)))
  } #end of MatchBalance


balanceMV  <- function(formul, data=NULL, match.out=NULL, maxit=1000, weights=rep(1,nrow(data)),
                       nboots=100, nmc=nboots, verbose=0, ...)
  {

    tol <- .Machine$double.eps*100
    
    if(!is.list(match.out) & (!is.null(match.out))) {
      warning("'Match' object contains no valid matches")
      match.out  <- NULL
    }

    if ( (class(match.out) != "Match") & (!is.null(match.out)) ) {
      warning("Object not of class 'Match'")
      match.out  <- NULL
    }

    if (is.null(data))
      {
        data  <- as.data.frame(model.frame(formul))
        xdata <- get.xdata(formul,datafr=environment(formul))
        Tr    <- get.ydata(formul,datafr=environment(formul))
      } else {
        data  <- as.data.frame(data)
        xdata <- get.xdata(formul, data)
        Tr <- get.ydata(formul, data)
      }

    if (sum(is.na(data)!=0))
      stop("MatchBalance: NAs found in data input")

    if (nboots!=0 & nboots < 10)
      {
        nboots <- 10
        warning("at least 10 nboots must be run")
      }

    if ( (nmc!=0 & nmc < 10) | (nboots > 0 & nmc < 10) )
      {
        nmc <- 10
        warning("at least 10 Monte Carlo runs are required")
      }

    MC <- FALSE
    if (nmc > 0)
      MC <- TRUE
    
    if (!is.null(data$weight37172))
      stop("balanceMV(): data$weight37172 cannot exist")

    data$weight37172  <- weights

    ddata  <- list()
    ddata$Tr <- Tr
    ddata$nboots <- nboots
    ddata$weights.unmatched <- data$weight37172        

    #FULL SAMPLE RUN
    t1  <- glm(formul, family=binomial(link=logit), weights=weight37172,
               control=glm.control(maxit=maxit, ...), data=data)

    t1$df <- t1$df.null-t1$df.residual
    t1$stat <-  (t1$null.deviance-t1$deviance)
    t1$p.value  <- 1-pchisq((t1$null.deviance-t1$deviance),t1$df)

    t1$treated.fitted  <- t1$fitted.values[Tr==1]
    t1$control.fitted  <- t1$fitted.values[Tr==0]

    bm.ks  <- Mks.test(t1$treated.fitted, t1$control.fitted, MC=MC)

    #Monte Carlo Estimates for the full sample before matching ks.test
    Ys <- t1$fitted
    obs <- nrow(data)
    cutp <- round(obs/2)
    sim.ks.bm.pval <- NULL
    bbcount <- 0
    if (nmc > 0)
      {
        for (bb in 1:nmc)
          {
            sindx  <- sample(1:obs, obs, replace=TRUE)
            
            
            X1tmp  <- Ys[sindx[1:cutp]]
            X2tmp  <- Ys[sindx[(cutp+1):obs]]
            
            s.ks   <- Mks.test(X1tmp, X2tmp, exact=FALSE, MC=MC)$statistic
            
            if (s.ks >= (bm.ks$statistic - tol) )
              bbcount  <- bbcount + 1
            
          }
        sim.ks.bm.pval  <- bbcount/nmc
      }

    #BOOTSTRAP SAMPLE RUN
    if (nboots > 0)
      {
        bcount <- 0
        if (verbose > 1)
          cat("unmatched bootstraps:\n")

        na.indx   <- is.na(t1$coeff)
        use.coeff <- t1$coeff[!na.indx]
        Xs        <- xdata[,!na.indx]            
        sigma     <- summary(t1)$cov.scaled
        
        parms.s <- MASS::mvrnorm(nboots, mu=use.coeff, Sigma=sigma)
          
        for (s in 1:nboots)
          {
            if (verbose > 1)
              cat("s: ",s,"\n")

            mu <- Xs %*% as.matrix(parms.s[s,])
            Ys <- exp(mu)/(1 + exp(mu))

            bbcount <- 0
            for (bb in 1:nmc)
              {
                sindx  <- sample(1:obs, obs, replace=TRUE)

                X1tmp  <- Ys[sindx[1:cutp]]
                X2tmp  <- Ys[sindx[(cutp+1):obs]]

                s.ks   <- Mks.test(X1tmp, X2tmp, exact=FALSE, MC=MC)

                if (s.ks$statistic >= (bm.ks$statistic - tol) )
                  bbcount  <- bbcount + 1                
              } #end of nmc
            bb.pval1  <- bbcount/nmc

            if (bb.pval1 >= sim.ks.bm.pval )
              {
                bcount  <- bcount+1
              }

            if ((bb.pval1==sim.ks.bm.pval) & (sim.ks.bm.pval==0))
              {
                bcount  <- bcount-1
              }
            
          }#end of s loop
        p.bm.boot <- bcount/nboots        
      } else {
        p.bm.boot <- NULL
      } #end of if nboots


    if(!is.null(match.out))
      {
        Mdata  <- rbind(as.data.frame(data[match.out$index.treated,]), data[match.out$index.control,])
        Mdata$weight37172  <- c(match.out$weights*weights[match.out$index.treated], 
                                match.out$weights*weights[match.out$index.control])
	MTr <- c(rep(1,length(match.out$index.treated)), rep(0,length(match.out$index.control)))

        xdata <- as.data.frame(get.xdata(formul, Mdata))
        xdata <- as.matrix(xdata)
        Xs    <- xdata
        intercept <- FALSE
        if (sum(xdata[,1]==rep(1,nrow(xdata)))==nrow(xdata))
          {
            intercept <- TRUE
            xdata <- xdata[,2:ncol(xdata)]

            Mt1  <- glm(MTr~xdata, family=quasibinomial(link=logit), 
                        control=glm.control(maxit=maxit, ...),
                        weights=Mdata$weight37172)
          } else {
            Mt1  <- glm(MTr~xdata-1, family=quasibinomial(link=logit), weights=Mdata$weight37172,
                        control=glm.control(maxit=maxit, ...))                        
          }

        Mdata0  <- list()
        Mdata0$y  <- MTr
        Mdata0$weight37172  <- Mdata$weight37172
        
        Mt0  <- glm(y~1, family=quasibinomial(link=logit), weights=Mdata0$weight37172,
                    control=glm.control(maxit=maxit, ...), data=Mdata0)
        
        rd0<-residuals(Mt0) # deviance residuals
        rd1<-residuals(Mt1) # deviance residuals
        Mt1$df <- Mt1$df.null-Mt1$df.residual
        #the square root of the weights is already included in the residuals
        wdeviance1  <- sum( rd1^2 *sqrt(Mdata$weight37172) )
        wdeviance0  <- sum( rd0^2 *sqrt(Mdata$weight37172) )
        rr  <- (rd1^2)*sqrt(Mdata$weight37172)    
        Mt1$dispersion <- sum( rr )/( sum(Mdata$weight37172) - Mt1$df)
        Mt1$stat  <- (wdeviance0-wdeviance1)/Mt1$dispersion
        Mt1$p.value  <- 1-pchisq(Mt1$stat,Mt1$df)

        Mt1$treated.fitted  <- Mt1$fitted.values[MTr==1]
        Mt1$control.fitted  <- Mt1$fitted.values[MTr==0]
        
        am.ks  <- Mks.test(Mt1$treated.fitted, Mt1$control.fitted, MC=MC)

	ddata$Tr.matched <- MTr        
	ddata$weights.matched <- Mdata$weight37172
	ddata$weights <- match.out$weights
        ddata$index.treated.indata  <- match.out$index.treated.indata
        ddata$index.control.indata  <- match.out$index.control.indata

        #Monte Carlo Estimates for the full sample before matching ks.test
        Ys <- Mt1$fitted
        obs <- nrow(Mdata)
        cutp <- round(obs/2)
        sim.ks.am.pval <- NULL
        bbcount <- 0
        if (nmc > 0)
          {
            for (bb in 1:nmc)
              {
                sindx  <- sample(1:obs, obs, replace=TRUE)
                
                X1tmp  <- Ys[sindx[1:cutp]]
                X2tmp  <- Ys[sindx[(cutp+1):obs]]
                
                s.ks   <- Mks.test(X1tmp, X2tmp, exact=FALSE, MC=MC)$statistic
                
                if (s.ks >= (am.ks$statistic - tol) )
                  bbcount  <- bbcount + 1
                
              }
            sim.ks.am.pval  <- bbcount/nmc
          }

        #BOOTSTRAP SAMPLE RUN
        if (nboots > 0)
          {
            bcount <- 0
            
            if (verbose > 1)
              cat("matched bootstraps:\n")

            na.indx   <- is.na(Mt1$coeff)
            use.coeff <- Mt1$coeff[!na.indx]
            Xs        <- Xs[,!na.indx]            
            sigma     <- summary(Mt1, dispersion=Mt1$dispersion)$cov.scaled
            
            parms.s <- MASS::mvrnorm(nboots, mu=use.coeff, Sigma=sigma)
          
            for (s in 1:nboots)
              {
                if (verbose > 1)
                  cat("s: ",s,"\n")
                
                mu <- Xs %*% as.matrix(parms.s[s,])
                Ys <- exp(mu)/(1 + exp(mu))

                bbcount <- 0
                for (bb in 1:nmc)
                  {
                    sindx  <- sample(1:obs, obs, replace=TRUE)
                    
                    X1tmp  <- Ys[sindx[1:cutp]]
                    X2tmp  <- Ys[sindx[(cutp+1):obs]]

                    s.ks   <- Mks.test(X1tmp, X2tmp, exact=FALSE, MC=MC)

                    if (s.ks$statistic >= (am.ks$statistic - tol) )
                      bbcount  <- bbcount + 1                
                  } #end of nmc
                bb.pval1  <- bbcount/nmc

                if (bb.pval1 >= sim.ks.am.pval)
                  {
                    bcount  <- bcount+1
                  }

                if ((bb.pval1==sim.ks.am.pval) & (sim.ks.am.pval==0))
                  {
                    bcount  <- bcount-1
                  }                
              }#end of s loop
            p.am.boot <- bcount/nboots        
          } else {
            p.am.boot <- NULL
          } #end of if nboots
        
      } else {
        am.ks  <- NULL
        Mt1  <- NULL
        p.am.boot <- NULL

      }

    z  <- list(pval.kboot.unmatched=p.bm.boot, pval.kboot.matched=p.am.boot,
               logit.unmatched=t1, logit.matched=Mt1,               
               ks.unmatched=bm.ks, ks.matched=am.ks,
               data=ddata)
    class(z)  <- "balanceMV"

    return(z)
  } #end of balanceMV

summary.balanceMV  <- function(object, ..., digits=5)
  {
    if (class(object) != "balanceMV") {
      warning("Object not of class 'balanceMV'")
      return(NULL)
    }

    t1.df <- object$logit.unmatched$df.null-object$logit.unmatched$df.residual
    t1.stat <-  object$logit.unmatched$null.deviance-object$logit.unmatched$deviance

    if(!is.null(object$ks.matched))
      {    
        Mt1.df  <- object$logit.matched$df.null-object$logit.matched$df.residual
        Mt1.stat  <- object$logit.matched$null.deviance-object$logit.matched$deviance
      }

    cat("\nBefore Matching:\n")
    w <- object$data$weights.unmatched[object$data$Tr==1]
    mean.tr.unmatched <- sum(object$logit.unmatched$fitted[object$data$Tr==1]*w)/sum(w)

    cat("Mean Probability of Treatment for Treated Observations:",
        format(mean.tr.unmatched,digits=digits),"\n")
    
    w <- object$data$weights.unmatched[object$data$Tr==0]
    mean.co.unmatched <- sum(object$logit.unmatched$fitted[object$data$Tr==0]*w)/sum(w)
    
    cat("Mean Probability of Treatment for Control Observations:",
        format(mean.co.unmatched,digits=digits),"\n")
    
    if(!is.null(object$ks.matched))
      {
        mean.tr.matched <- 
          sum(object$logit.matched$fitted[object$data$Tr.matched==1]*object$data$weights)/sum(object$data$weights)

        mean.co.matched <- 
          sum(object$logit.matched$fitted[object$data$Tr.matched==0]*object$data$weights)/sum(object$data$weights)

        cat("\nAfter Matching:\n")    
        cat("Mean Probability of Treatment for Treated Observations:", format(mean.tr.matched, digits=digits),"\n")
        cat("Mean Probability of Treatment for Control Observations:", format(mean.co.matched, digits=digits),"\n")
      }

    cat("\nKolmogorov-Smirnov Test for Balance Before Matching:\n")
    cat("statistic:",format(object$ks.unmatched$statistic,digits=digits), " p-val: ",
        format.pval(object$ks.unmatched$p.value,digits=digits),"\n")
    if( object$data$nboots > 0)
      cat("bootstrap p-value:",format.pval(object$pval.kboot.unmatched,digits=digits),"\n")      
      

    if(!is.null(object$ks.matched))
      {    
        cat("\nKolmogorov-Smirnov Test for Balance After Matching:\n")
        cat("statistic:",format(object$ks.matched$statistic,digits=digits), " p-val: ",
            format.pval(object$ks.matched$p.value,digits=digits),"\n")
        if( object$data$nboots > 0)
          cat("bootstrap p-value:",format.pval(object$pval.kboot.matched,digits=digits),"\n")      
      }

    cat("\nChi-Square Deviance Test for Balance Before Matching:\n")
    cat("statistic:",format(t1.stat,digits=digits),
        " p-val:",format.pval(object$logit.unmatched$p.value,digits=digits),
        "df:",format(t1.df,digits=digits),"\n")

    if(!is.null(object$ks.matched))
      {
        cat("\nChi-Square Deviance Test for Balance After Matching:\n")
        cat("statistic:",format(Mt1.stat,digits=digits),
            " p-val:",format.pval(object$logit.matched$p.value,digits=digits),"df:",
            format(Mt1.df,digits=digits),"\n")
      }

    cat("\n\n")
  }#end of summary.Match


get.xdata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "term.labels"))==0 & attr(t1, "intercept")==0) {
    m <- NULL;  # no regressors specified for the model matrix
  } else {
    m <- model.matrix(formul, data=datafr, drop.unused.levels = TRUE);
  }
  return(m);
}


# get.ydata:
# Return response vector corresponding to the formula in formul
# 
get.ydata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "response"))==0) {
    m <- NULL;  # no response variable specified
  }  else {
    m <- model.response(model.frame(formul, data=datafr));
  }
  return(m);
}

#
# bootstrap ks test implemented
#

ks.boot  <- function(Tr, Co, nboots=1000, verbose=0)
  {
    tol <- .Machine$double.eps*100
    Tr <- Tr[!is.na(Tr)]
    Co <- Co[!is.na(Co)]

    w    <- c(Tr, Co)
    obs  <- length(w)
    cutp <- round(obs/2)
    ks.boot.pval <- NULL
    bbcount <- 0

    if (nboots < 10)
      {
        nboots  <- 10
        warning("At least 10 'nboots' must be run; seting 'nboots' to 10")
      }

    if (nboots < 500)
      warning("For publication quality p-values it is recommended that 'nboots'\n be set equal to at least 500 (preferably 1000)") 
    
    fs.ks  <- Mks.test(Tr, Co, MC=TRUE)    

    for (bb in 1:nboots)
      {
        if (verbose > 1)
          cat("s:", bb, "\n")
        
        sindx  <- sample(1:obs, obs, replace=TRUE)
        
        
        X1tmp  <- w[sindx[1:cutp]]
        X2tmp  <- w[sindx[(cutp+1):obs]]
        
        s.ks   <- Mks.test(X1tmp, X2tmp, exact=FALSE, MC=TRUE)$statistic
        
        if (s.ks >= (fs.ks$statistic - tol) )
          bbcount  <- bbcount + 1
        
      }
    ks.boot.pval  <- bbcount/nboots

    ret  <- list(ks.boot.pvalue=ks.boot.pval, ks=fs.ks, nboots=nboots)
    class(ret)  <- "ks.boot"

    return(ret)
  } #end of ks.boot

summary.ks.boot <- function(object, ..., digits=5)
  {
    if(!is.list(object)) {
      warning("object not a valid 'ks.boot' object")
      return()
    }

    if (class(object) != "ks.boot") {
      warning("Object not of class 'ks.boot'")
      return()
    }    
        
    cat("\n")
    cat("Bootstrap p-value:    ", format.pval(object$ks.boot.pvalue, digits=digits), "\n")
    cat("Full Sample Statistic:", format(object$ks$statistic, digits=digits), "\n")
    cat("Naive p-value:        ", format(object$ks$p.value, digits=digits), "\n")
#    cat("nboots completed      ", object$nboots, "\n")
#   cat("\n")
  } #end of summary.ks.boot




