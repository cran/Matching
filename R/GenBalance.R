Mt.test.pvalue  <- function(Tr, Co, weights)
  {
    v1  <- Tr-Co
    estimate  <- sum(v1*weights)/sum(weights)
    var1  <- sum( ((v1-estimate)^2)*weights )/( sum(weights)*sum(weights) )
    
    statistic  <- estimate/sqrt(var1)
#    p.value    <- (1-pnorm(abs(statistic)))*2
    p.value    <- (1-pt(abs(statistic), df=sum(weights)-1))*2

    return(p.value)
  } #end of Mt.test.pvalue

GenBalance <- function(rr, X, 
                       nvars=ncol(X), nboots = 0, ks=TRUE, verbose = FALSE)
{
  
  index.treated <- rr[,1]
  index.control <- rr[,2]
  weights <- rr[,3]

  tol  <- .Machine$double.eps*100  
  storage.t <- c(rep(NA,nvars))
  storage.k <- c(rep(NA,nvars))
  fs.ks     <- matrix(nrow=nvars, ncol=1)
  s.ks      <- matrix(nrow=nvars, ncol=1)  
  bbcount   <- matrix(0, nrow=nvars, ncol=1)
  dummy.indx  <- matrix(0, nrow=nvars, ncol=1)
  
  w  <- c(X[,1][index.treated], X[,1][index.control])
  obs <- length(w)
  n.x  <- length(X[,1][index.treated])
  n.y  <- length(X[,1][index.control])
  cutp <- round(obs/2)  
  w  <- matrix(nrow=obs, ncol=nvars)

  for (i in 1:nvars)
    {
      w[,i] <- c(X[,i][index.treated], X[,i][index.control])
      
      t.out <- Mt.test.pvalue(X[,i][index.treated],
                              X[,i][index.control],
                              weights = weights)
      
      storage.t[i] <- t.out            
      
#      print(length(unique(X[,i])) < 3)

      
      dummy.indx[i]  <- length(unique(X[,i])) < 3

      if (!dummy.indx[i] & ks & nboots > 9)
        {
          fs.ks[i]  <- ks.fast(X[,i][index.treated], X[,i][index.control],
                               n.x=n.x, n.y=n.y, n=obs)
        } else if(!dummy.indx[i] & ks)
          {
            #storage.k[i] <- ks.test(X[,i][index.treated], X[,i][index.control])$p.value
            
            storage.k[i] <- Mks.test(X[,i][index.treated], X[,i][index.control],
                                     MC=TRUE)$p.value

          }
    }#end of i loop

  if (ks & nboots > 9)
    {
      n.x  <- cutp
      n.y  <- obs-cutp
      for (b in 1:nboots)
        {
          sindx  <- sample(1:obs, obs, replace = TRUE)
          
          for (i in 1:nvars)
            {
              
              if (dummy.indx[i])
                next;
              
              X1tmp <- w[sindx[1:cutp],i ]
              X2tmp <- w[sindx[(cutp + 1):obs], i]
              s.ks[i] <- ks.fast(X1tmp, X2tmp, n.x=n.x, n.y=n.y, n=obs)
              if (s.ks[i] >= (fs.ks[i] - tol) )
                bbcount[i]  <-  bbcount[i] + 1
            }#end of i loop
        } #end of b loop
      
      for (i in 1:nvars)
        {
          
          if (dummy.indx[i])
            {
              storage.k[i]  <- 999
              next;
            }
          
          storage.k[i]  <- bbcount[i]/nboots
          
        }

#  cat("bbcount:\n")
#  print(bbcount)
#  cat("storage:\n")
#  print(storage.k)

      output <- min(c(storage.t, storage.k), na.rm = TRUE)
    } else if(ks){
      output <- min(c(storage.t, storage.k), na.rm = TRUE)
    } else {
      output <- min(storage.t, na.rm = TRUE)
    }
  
  if (verbose == TRUE)
    {
      for (i in 1:nvars)
        {
          cat("\n", i, " t-test p-value =", storage.t[i], "\n" )
          cat("\n", i, " boot ks-test p-value =", storage.k[i], "\n" )
        }
      cat("\nreturn value:", output, "\n\n")
    }
  
  return(output)
} #end of genBalance

#
# writing fast KS test
#

ks.fast <- function(x, y, n.x, n.y, n)
  {
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]        
    
    return( max(abs(z)) )
  }


KSbootBalanceSummary <- function(index.treated, index.control, X, 
                                 nboots = 1000)
{
  nvars <- ncol(X)

  tol  <- .Machine$double.eps*100  
  storage.k <- c(rep(NA,nvars))
  storage.k.naive <- c(rep(NA,nvars))  
  fs.ks     <- matrix(nrow=nvars, ncol=1)
  s.ks      <- matrix(nrow=nvars, ncol=1)  
  bbcount   <- matrix(0, nrow=nvars, ncol=1)
  dummy.indx  <- matrix(0, nrow=nvars, ncol=1)
  
  w  <- c(X[,1][index.treated], X[,1][index.control])
  obs <- length(w)
  n.x  <- length(X[,1][index.treated])
  n.y  <- length(X[,1][index.control])
  cutp <- round(obs/2)  
  w  <- matrix(nrow=obs, ncol=nvars)

  for (i in 1:nvars)
    {
      w[,i] <- c(X[,i][index.treated], X[,i][index.control])
      
      dummy.indx[i]  <- length(unique(X[,i])) < 3

      if (!dummy.indx[i])
        {
          foo <- Mks.test(X[,i][index.treated], X[,i][index.control], MC=TRUE)
          fs.ks[i] <- foo$statistic
          storage.k.naive[i] <- foo$p.value
        } 
    }#end of i loop

  n.x  <- cutp
  n.y  <- obs-cutp
  for (b in 1:nboots)
    {
      sindx  <- sample(1:obs, obs, replace = TRUE)
      
      for (i in 1:nvars)
        {
          
          if (dummy.indx[i])
            next;
          
          X1tmp <- w[sindx[1:cutp],i ]
          X2tmp <- w[sindx[(cutp + 1):obs], i]
          s.ks[i] <- ks.fast(X1tmp, X2tmp, n.x=n.x, n.y=n.y, n=obs)
          if (s.ks[i] >= (fs.ks[i] - tol) )
            bbcount[i]  <-  bbcount[i] + 1
        }#end of i loop
    } #end of b loop

  for (i in 1:nvars)
    {
      if (dummy.indx[i])
        {
          storage.k[i]  <- NA
          next;
        }
      
      storage.k[i]  <- bbcount[i]/nboots
    }

  ret = list(ks.boot.pval=storage.k, ks.naive.pval=storage.k.naive, ks.stat=fs.ks)

  return(ret)
} #end of KSbootBalanceSummary
