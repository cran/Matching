Matchby  <- function(Y, Tr, X, by, estimand="ATT", M=1, ties=FALSE, replace=TRUE,
                     exact = NULL, caliper = NULL, 
                     Weight = 1, Weight.matrix = NULL,
                     tolerance = 1e-05, distance.tolerance = 1e-05,                    
                     print.level=1, version="fast", ...)
  {

    X  <- as.matrix(X)
    xvars  <- ncol(X)

    #index for raw obs
    nobs  <- length(Tr)
    nl  <- 1:nobs

    orig.treated.nobs  <- sum(Tr==1)
    t.dta  <- cbind(Y, Tr, X)
    ncols  <- ncol(t.dta)
    t.dta  <- split(t.dta, by, drop=TRUE)

    nindx  <- length(t.dta)

    t.treated  <- NULL
    t.control  <- NULL
    t.w  <- NULL

    if (replace!=FALSE & replace!=TRUE)
      {
        warning("'replace' must be TRUE or FALSE.  Setting to TRUE")
        replace <- TRUE
      }
    if(replace==FALSE)
      {
        ties <- FALSE        
      }    
    if (ties!=FALSE & ties!=TRUE)
      {
        warning("'ties' must be TRUE or FALSE.  Setting to TRUE")
        ties <- TRUE
      }    
    
    ndrops  <- 0
    for (i in 1:nindx)
      {
        if(print.level > 0)
          cat(i,"of", nindx, "groups\n");
        
        f.dta  <- matrix(t.dta[[i]], ncol=ncols)

        if(nrow(f.dta) < 2)
          {
            ndrops  <- ndrops + nrow(f.dta)
            next;
          }
        
        if( var(f.dta[,2])==0 )
          {
            ndrops  <- ndrops + nrow(f.dta)            
            next;
          }

        X = f.dta[,3:ncols]
        if (xvars > 1)
          {
            indx  <- apply(X, 2, var) < tolerance

            if (sum(indx) == ncols)
              {
                ndrops  <- ndrops + nrow(f.dta)                
                next;
              }

            X = X[,!indx]
          } else {
            if (var(X) < tolerance)
              {
                ndrops  <- ndrops + nrow(f.dta)                
                next;
              }
          }

        t1  <- Match(Y=f.dta[,1], Tr=f.dta[,2], X=X,
                     estimand=estimand, M=M,
                     exact=exact, caliper=caliper, replace=replace, ties=ties,
                     Weight=Weight, Weight.matrix=Weight.matrix,
                     tolerance=tolerance, distance.tolerance=distance.tolerance,
                     version=version, ...)

        ndrops <- ndrops + t1$ndrops.matches 

        t.treated  <- c(t.treated, t1$mdata$Y[t1$mdata$Tr==1])
        t.control  <- c(t.control, t1$mdata$Y[t1$mdata$Tr==0])
        t.w  <- c(t.w, t1$weights)
      }
                   
    sum.tw  <- sum(t.w)
    est  <- sum(t.treated*t.w)/sum.tw  - sum(t.control*t.w)/sum.tw
    v1  <- t.treated - t.control
    varest  <- sum( ((v1-est)^2)*t.w)/(sum.tw*sum.tw)
    se.standard  <- sqrt(varest)
    
    ret  <- list(est=est, se.standard=se.standard,
                          ret=cbind(t.treated, t.control, t.w))
    ret$orig.nobs  <- nobs
    ret$nobs  <- length(t.treated)
    ret$wnobs  <- sum.tw
    ret$orig.treated.nobs  <- orig.treated.nobs
    ret$ndrops  <- ndrops
    ret$estimand  <- estimand
    ret$version  <- version
    class(ret)  <- "Matchby"
    
    invisible(return(ret))
  }

summary.Matchby  <- function(object, ..., digits=5)
  {
    if(!is.list(object)) {
      warning("'Matchby' object contains less than two valid matches.  Cannot proceed.")
      return(invisible(NULL))
    }
    
    if (class(object) != "Matchby") {
      warning("Object not of class 'Matchby'")
      return(invisible(NULL))
    }

    cat("\n")
    cat("Estimate... ",format(object$est,digits=digits),"\n")
    cat("SE......... ",format(object$se.standard,digits=digits),"\n")
    cat("T-stat..... ",format(object$est/object$se.standard,digits=digits),"\n")
    cat("p.val...... ",format.pval((1-pnorm(abs(object$est/object$se.standard)))*2,digits=digits),"\n")
    cat("\n")        

    cat("Original number of observations.............. ", object$orig.nobs,"\n")
    if(object$estimand!="ATC")
      {
        cat("Original number of treated obs............... ", object$orig.treated.nobs,"\n")
      } else {
        cat("Original number of control obs............... ",
            object$orig.nobs- object$orig.treated.nobs,"\n")        
      }
    cat("Matched number of observations............... ", round(object$wnobs, 3),"\n")
    cat("Matched number of observations  (unweighted). ", object$nobs,"\n")
    cat("\n")
        
    if ( object$ndrops > 0)
        {
          cat("Number of observations dropped............... ", object$ndrops, "\n")
          cat("\n")
        }
  } #end of summary.Matchby

