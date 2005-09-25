/* #include <sys/time.h> */

#include "scythematrix.h"

using namespace SCYTHE;
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

#include "matching.h"

extern "C"
{
#include <Rdefines.h>

  SEXP FastMatchC(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
		  SEXP I_ww, SEXP I_Tr, SEXP I_X, SEXP I_weight)
  {
    SEXP ret;
    
    long N, xvars, All, M;
    double cdd;

    long i, j, k, r, c;

    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);

    //    struct timeval tv1, tv2;	// Required "timeval" structures (see man pg.)
    //    struct timezone tz1, tz2;	// Required "timezone" structures (see man pg.)
    //    long tret = gettimeofday(&tv1,&tz1);

    Matrix ww = Matrix(xvars, xvars);
    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    Matrix weight = Matrix(N, 1);

    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<xvars; j++)
	{
	for(i=0;i<xvars; i++)
	  {
	    //ww(i,j) = REAL(I_ww)[k];
	    ww[M(i,j,xvars)] = REAL(I_ww)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	Tr[i] = REAL(I_Tr)[i];
      }

    //rows and colums are fliped!! j,i != i,j
    k=0;
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<N; i++)
	  {
	    //X(i,j) = REAL(I_X)[k];
	    X[M(i,j,xvars)] = REAL(I_X)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	weight[i] = REAL(I_weight)[i];
      }

    Matrix INN = seqa(1, 1, N);
    // TT is just equal to INN
    Matrix index_onesN = ones(N, 1);

    Matrix Kcount  = zeros(N,1);
    Matrix KKcount = zeros(N,1);
    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0, ACTMATsum = 0;

    Matrix xx;
    Matrix DX;
    Matrix Dist, DistPot, weightPot, tt, 
      weightPot_sort, weightPot_sum, ACTDIST, Wi, xmat, I, IM, W;
    Matrix ACTMAT, POTMAT;

    //    tret = gettimeofday(&tv2,&tz2);
    //    long secs = tv2.tv_sec - tv1.tv_sec;
    //    long msecs = tv2.tv_usec - tv1.tv_usec;
    //    double actual = ((double) secs*1000000+ (double) msecs)/1000000;
    //    printf("actual: %lf, secs: %d, msecs: %d\n", actual, secs, msecs);

    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    double *S = (double *) malloc(N*sizeof(double));  

    // struct timeval tv1, tv2;	// Required "timeval" structures (see man pg.)
    // struct timezone tz1, tz2;	// Required "timezone" structures (see man pg.)
    // long tret = gettimeofday(&tv1,&tz1);

    int first=1;
    for(i=0; i < N; i++)
      {
	// treatment indicator for observation to be matched        
        TREATi = (int) Tr[i];
	
        // proceed with all observations if All==1
        // but only with treated observations if All=0        
        if ( (TREATi==1 & All!=1) | All==1 )
          {
            // covariate value for observation to be matched                        
            xx = X(i,_);


            DX = (X - (index_onesN * xx)) * t(ww);

            if (xvars>1)
              {
                //JSS
                foo1 = t(multi_scalar(DX, DX));
		Dist = t(sumc(foo1));

	      } 
	    else 
	      {
		Dist = multi_scalar(DX, DX);
	      } // end of xvars

	    // Dist distance to observation to be matched
            // is N by 1 vector	    

            // set of potential matches (all observations with other treatment)
            // JSS, note:logical vector
	    POTMAT = EqualityTestScalar(Tr, 1-TREATi);

            // X's for potential matches
            DistPot = selif(Dist, POTMAT);
            weightPot = selif(weight, POTMAT);

	    long weightPot_size = size(weightPot);

	    for(j=0; j< weightPot_size; j++)
	      {
		// assume that weightPot_size = size(DistPot)
		order_DistPot[j] = j;
		S[j] = (double) DistPot[j];
	      }

	    rsort_with_index (S, order_DistPot, weightPot_size);

	    weightPot_sort = Matrix(weightPot_size, 1);
	    for(j=0; j < weightPot_size; j++)
	      {
		weightPot_sort[j] = weightPot[order_DistPot[j]];
	      }
            weightPot_sum = cumsum(weightPot_sort);

	    tt = seqa(1, 1, rows(weightPot_sum));

	    foo1 = GreaterEqualTestScalar(weightPot_sum, M);
	    foo2 = selif(tt, foo1);

	    long MMM = (long) min(foo2) - 1;

	    // distance at last match
            double Distmax = S[MMM];

            // selection of actual matches 
            // logical index
	    ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	    ACTMAT = VectorAnd(POTMAT, ACTMAT);

            // distance to actual matches            
            ACTDIST = selif(Dist, ACTMAT);

            // counts how many times each observation is matched.
	    double Mi = sum(multi_scalar(weight, ACTMAT));

            Kcount = Kcount + 
	      weight[i] * multi_scalar(weight, ACTMAT)/Mi;

	    foo1 = multi_scalar(weight, weight);
	    foo1 = multi_scalar(foo1, ACTMAT);

	    KKcount = KKcount + (weight[i]* foo1)/(Mi*Mi);

            Wi = selif(weight, ACTMAT);
	    Wi = weight[i]*Wi/Mi;

	    ACTMATsum = (int) sumc(ACTMAT)[0];

	    // collect results
	    if (first==1)
	      {
		I = ones(ACTMATsum, 1)*(i+1);
		IM = selif(INN, ACTMAT);
		W = Wi;
		first = 0;
	      }// end of first==1 
	    else 
	      {
		I = rbind(I, ones(ACTMATsum, 1)*(i+1));
		IM = rbind(IM, selif(INN, ACTMAT));
		W = rbind(W, Wi);
	      } // end of i else

	  } // end of (TREATi==1 & All!=1) | All==1 )
      } //END OF i MASTER LOOP!

    // tret = gettimeofday(&tv2,&tz2);
    // long secs = tv2.tv_sec - tv1.tv_sec;
    // long msecs = tv2.tv_usec - tv1.tv_usec;
    // double actual = ((double) secs*1000000+ (double) msecs)/1000000;
    // printf("actual: %lf, secs: %d, msecs: %d\n", actual, secs, msecs);

    /*ATT is okay already */
    /* ATE*/
    if(All==1)
      {
	long tl  = rows(I);
	Matrix I2  = zeros(tl, 1);
	Matrix IM2 = zeros(tl, 1);
	Matrix trt = zeros(tl, 1);

	for(i=0; i<tl; i++)
	  {
	    k =(int) I[i] -1 ;
	    trt[i] = Tr[k];
	  }

	for (i=0; i<tl; i++)
	  {
	    if (trt[i]==1)
	      {
		I2[i] = I[i];
		IM2[i] = IM[i];
	      } 
	    else
	      {
		I2[i] = IM[i];
		IM2[i] = I[i];		
	      }
	  }

	I = I2;
	IM = IM2;
      } 
    else if(All==2)     /* ATC */
      {
	Matrix Itmp = I;
	Matrix IMtmp = IM;

	I = IMtmp;
	IM = Itmp;
      }

    Matrix rr = cbind(I, IM);
    rr = cbind(rr, W);

    // display(rr);

    /* Free Memory */
    free(S);
    free(order_DistPot);

    r = rows(rr);
    c = cols(rr);
    
    PROTECT(ret=allocMatrix(REALSXP, r, c));
    /* Loop through the data and display the same in matrix format */
    k = 0;
    for( i = 0; i < c; i++ )
      {
	for( j = 0; j < r; j++ )
	  {
	    // REAL(ret)[k] = rr(j,i);
	    // REAL(ret)[k] = rr[j*c+i];
	    /* Use Macro to Index */
	    REAL(ret)[k] = rr[M(j, i, c)];
	    k++;
	  }
      }
    UNPROTECT(1);
    return(ret);    
  } //end of FastMatchC

  SEXP MatchLoopC(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
		  SEXP I_caliper,
		  SEXP I_ww, SEXP I_Tr, SEXP I_X, SEXP I_weight,
		  SEXP I_CaliperVec, SEXP I_Xorig,
		  SEXP I_restrict_trigger, SEXP I_restrict_nrow, SEXP I_restrict)
  {
    SEXP ret;

    long N, xvars, All, M, caliper, restrict_trigger, restrict_nrow, sum_caliper_drops=0;
    double cdd, diff, foo_d;

    long i, j, k, r, c;

    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    caliper = (long) asReal(I_caliper);
    restrict_nrow = asInteger(I_restrict_nrow);
    restrict_trigger = asInteger(I_restrict_trigger);

    //    struct timeval tv1, tv2;	// Required "timeval" structures (see man pg.)
    //    struct timezone tz1, tz2;	// Required "timezone" structures (see man pg.)
    //    long tret = gettimeofday(&tv1,&tz1);

    Matrix ww = Matrix(xvars, xvars);
    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    Matrix weight = Matrix(N, 1);

    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<xvars; j++)
	{
	for(i=0;i<xvars; i++)
	  {
	    //ww(i,j) = REAL(I_ww)[k];
	    ww[M(i,j,xvars)] = REAL(I_ww)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	Tr[i] = REAL(I_Tr)[i];
      }

    //rows and colums are fliped!! j,i != i,j
    k=0;
    for(j=0;j<xvars; j++)
      {
	for(i=0;i<N; i++)
	  {
	    //X(i,j) = REAL(I_X)[k];
	    X[M(i,j,xvars)] = REAL(I_X)[k];
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	weight[i] = REAL(I_weight)[i];
      }

    Matrix IMi;
    Matrix Xorig;
    Matrix CaliperVec;
    if(caliper==1)
      {
	Xorig = Matrix(N, xvars);
	CaliperVec = Matrix(xvars, 1);

	for (i=0; i<xvars; i++)
	  {
	    CaliperVec[i] = REAL(I_CaliperVec)[i];
	  }

	//rows and colums are fliped!! j,i != i,j
	k=0;
	for(j=0;j<xvars; j++)
	  {
	    for(i=0;i<N; i++)
	      {
		//X(i,j) = REAL(I_X)[k];
		Xorig[M(i,j,xvars)] = REAL(I_Xorig)[k];
		k++;
	      }
	  }	
      } // end of caliper==1

    Matrix restrict; 
    if (restrict_trigger==1)
      {
	restrict = Matrix(restrict_nrow, 3);

	k=0;
	for(j=0;j<3; j++)
	  {
	    for(i=0;i<restrict_nrow; i++)
	      {
		restrict[M(i,j,3)] = REAL(I_restrict)[k];
		k++;
	      }
	  }	
      } /* if (restrict_trigger==1) */

    Matrix INN = seqa(1, 1, N);
    // TT is just equal to INN
    Matrix index_onesN = ones(N, 1);

    Matrix Kcount  = zeros(N,1);
    Matrix KKcount = zeros(N,1);
    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0, ACTMATsum = 0;

    Matrix xx;
    Matrix DX;
    Matrix Dist, DistPot, weightPot, tt, 
      weightPot_sort, weightPot_sum, ACTDIST, Wi, xmat, I, IM, W;
    Matrix ACTMAT, POTMAT;

    //    tret = gettimeofday(&tv2,&tz2);
    //    long secs = tv2.tv_sec - tv1.tv_sec;
    //    long msecs = tv2.tv_usec - tv1.tv_usec;
    //    double actual = ((double) secs*1000000+ (double) msecs)/1000000;
    //    printf("actual: %lf, secs: %d, msecs: %d\n", actual, secs, msecs);

    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    double *S = (double *) malloc(N*sizeof(double));  

    // struct timeval tv1, tv2;	// Required "timeval" structures (see man pg.)
    // struct timezone tz1, tz2;	// Required "timezone" structures (see man pg.)
    // long tret = gettimeofday(&tv1,&tz1);

    int first=1;
    for(i=0; i < N; i++)
      {
	// treatment indicator for observation to be matched        
        TREATi = (int) Tr[i];
	
        // proceed with all observations if All==1
        // but only with treated observations if All=0        
        if ( (TREATi==1 & All!=1) | All==1 )
          {
            // covariate value for observation to be matched                        
            xx = X(i,_);


            DX = (X - (index_onesN * xx)) * t(ww);

            if (xvars>1)
              {
                //JSS
                foo1 = t(multi_scalar(DX, DX));
		Dist = t(sumc(foo1));

	      } 
	    else 
	      {
		Dist = multi_scalar(DX, DX);
	      } // end of xvars

	    // Dist distance to observation to be matched
            // is N by 1 vector	    

	    if (restrict_trigger==1)
	      {
		for(j=0; j<restrict_nrow; j++)
		  {
		    if ( ((long) restrict[M(j,0,3)])-1 ==i )
		      {

                        if (restrict[M(j,2,3)] < 0) {
			  Dist[ ((long) restrict[M(j,1,3)])-1 ] = DOUBLE_XMAX;
			}
			else {
			    Dist[ ((long) restrict[M(j,1,3)])-1 ] = restrict[M(j,2,3)];
			}
		      }
		    else if ( ((long) restrict[M(j,1,3)])-1 ==i ) 
		      {

                        if (restrict[M(j,2,3)] < 0) {
			  Dist[ ((long) restrict[M(j,0,3)])-1 ] = DOUBLE_XMAX;
			}
			else {
                          Dist[ ((long) restrict[M(j,0,3)])-1 ] = restrict[M(j,2,3)];
			}
		      }
		  }
	      } /* if (restrict_trigger==1) */

            // set of potential matches (all observations with other treatment)
            // JSS, note:logical vector
	    POTMAT = EqualityTestScalar(Tr, 1-TREATi);

            // X's for potential matches
            DistPot = selif(Dist, POTMAT);
            weightPot = selif(weight, POTMAT);

	    long weightPot_size = size(weightPot);

	    for(j=0; j< weightPot_size; j++)
	      {
		// assume that weightPot_size = size(DistPot)
		order_DistPot[j] = j;
		S[j] = (double) DistPot[j];
	      }

	    rsort_with_index (S, order_DistPot, weightPot_size);

	    weightPot_sort = Matrix(weightPot_size, 1);
	    for(j=0; j < weightPot_size; j++)
	      {
		weightPot_sort[j] = weightPot[order_DistPot[j]];
	      }
            weightPot_sum = cumsum(weightPot_sort);

	    tt = seqa(1, 1, rows(weightPot_sum));

	    foo1 = GreaterEqualTestScalar(weightPot_sum, M);
	    foo2 = selif(tt, foo1);

	    long MMM = (long) min(foo2) - 1;

	    // distance at last match
            double Distmax = S[MMM];

            // selection of actual matches 
            // logical index
	    if (restrict_trigger==1)
	      {
		foo_d = min_scalar((Distmax+cdd), DOUBLE_XMAX-1);
		ACTMAT = LessEqualTestScalar(Dist,  foo_d);

		if ( ( (int) sumc(ACTMAT)[0]) < 1)
		  continue;
	      } 
	    else 
	      {
		ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	      } /* if (restrict_trigger==1) */

	    ACTMAT = VectorAnd(POTMAT, ACTMAT);

	    if (caliper==1)
	      {
		// counts how many times each observation is matched.
		long Mi = (long) sum(ACTMAT);

		IMi = selif(INN, ACTMAT);

		for (j=0; j<Mi; j++)
		  {
		    for (k=0; k<xvars; k++)
		      {
			diff = abs(Xorig[M(i, k, xvars)] - (double) Xorig[M((int) (IMi[j]-1),k,xvars)]); 
			if (diff > CaliperVec[k])
			  {
			    ACTMAT[(int) (IMi[j]-1)] = 0;
			    sum_caliper_drops++;
			    
			    break;
			  }
		      }
		  } 
		if ( ( (int) sumc(ACTMAT)[0]) < 1)
		  continue;
	      }//end of if caliper

            // distance to actual matches            
            ACTDIST = selif(Dist, ACTMAT);

            // counts how many times each observation is matched.
	    double Mi = sum(multi_scalar(weight, ACTMAT));

            Kcount = Kcount + 
	      weight[i] * multi_scalar(weight, ACTMAT)/Mi;

	    foo1 = multi_scalar(weight, weight);
	    foo1 = multi_scalar(foo1, ACTMAT);

	    KKcount = KKcount + (weight[i]* foo1)/(Mi*Mi);

            Wi = selif(weight, ACTMAT);
	    Wi = weight[i]*Wi/Mi;

	    ACTMATsum = (int) sumc(ACTMAT)[0];

	    // collect results
	    if (first==1)
	      {
		I = ones(ACTMATsum, 1)*(i+1);
		IM = selif(INN, ACTMAT);
		W = Wi;
		first = 0;
	      }// end of first==1 
	    else 
	      {
		I = rbind(I, ones(ACTMATsum, 1)*(i+1));
		IM = rbind(IM, selif(INN, ACTMAT));
		W = rbind(W, Wi);
	      } // end of i else

	  } // end of (TREATi==1 & All!=1) | All==1 )
      } //END OF i MASTER LOOP!

    // tret = gettimeofday(&tv2,&tz2);
    // long secs = tv2.tv_sec - tv1.tv_sec;
    // long msecs = tv2.tv_usec - tv1.tv_usec;
    // double actual = ((double) secs*1000000+ (double) msecs)/1000000;
    // printf("actual: %lf, secs: %d, msecs: %d\n", actual, secs, msecs);

    Matrix rr = cbind(I, IM);
    rr = cbind(rr, W);

    long tl  = rows(I);
    /*ATT is okay already */
    /* ATE*/
    if(All==1 & I[0]!=0)
      {
	Matrix I2  = zeros(tl, 1);
	Matrix IM2 = zeros(tl, 1);
	Matrix trt = zeros(tl, 1);

	for(i=0; i<tl; i++)
	  {
	    k =(int) I[i] -1 ;
	    trt[i] = Tr[k];
	  }

	for (i=0; i<tl; i++)
	  {
	    if (trt[i]==1)
	      {
		I2[i] = I[i];
		IM2[i] = IM[i];
	      } 
	    else
	      {
		I2[i] = IM[i];
		IM2[i] = I[i];		
	      }
	  }

	I = I2;
	IM = IM2;
      } 
    else if(All==2)     /* ATC */
      {
	Matrix Itmp = I;
	Matrix IMtmp = IM;

	I = IMtmp;
	IM = Itmp;
      }

    rr = cbind(rr, I);
    rr = cbind(rr, IM);

    if (caliper==1)
      {
	Matrix scalar_returns = zeros(tl, 1);
	scalar_returns[0] = sum_caliper_drops;
	rr = cbind(rr, scalar_returns);
      }

    // rr key
    // 1] I (unadjusted); 2] IM (unadjusted); 3] weight; 4] I (adjusted); 5] IM (adjusted);
    // scalar returns [0]: caliper drops

    /* Free Memory */
    free(S);
    free(order_DistPot);

    r = rows(rr);
    c = cols(rr);
    
    PROTECT(ret=allocMatrix(REALSXP, r, c));
    /* Loop through the data and display the same in matrix format */
    k = 0;
    for( i = 0; i < c; i++ )
      {
	for( j = 0; j < r; j++ )
	  {
	    // REAL(ret)[k] = rr(j,i);
	    // REAL(ret)[k] = rr[j*c+i];
	    /* Use Macro to Index */
	    REAL(ret)[k] = rr[M(j, i, c)];
	    k++;
	  }
      }
    UNPROTECT(1);
    return(ret);    
  } //end of MatchLoopC

} //end of extern "C"

double min_scalar (double a, double b)
{
  if (a < b)
    return(a);
  
  return(b);
} // end of min_scalar

double max_scalar (double a, double b)
{
  if (a > b)
    return(a);

  return(b);
} // end of max_scalar


Matrix multi_scalar (Matrix a, Matrix b)
{
  Matrix ret = a;
  long nrows = rows(a);
  long ncols = cols(a);

  for (long i = 0; i < nrows; i++) 
    {
      for (long j = 0; j < ncols; j++) 
	{
	  //ret(i, j) = a(i, j) * b(i, j);
	  ret[M(i, j, ncols)] = a[M(i, j, ncols)] * b[M(i, j, ncols)];
	}
    }

  return(ret);
} // multi_scalar


Matrix EqualityTestScalar(Matrix a, double s)
{
  long nrows = rows(a);
  long ncols = cols(a);
  Matrix ret =  zeros(nrows, ncols);

  for (long i =0; i< nrows; i++)
    {
      for (long j =0; j< ncols; j++)
	{
	  if( (a[M(i, j, ncols)] < (s+TOL)) &  (a[M(i, j, ncols)] > (s-TOL)) )
	    {
	      ret[M(i, j, ncols)] = 1;
	    }
	}
    }
  return(ret);
} //end of EqualityTestScalar

Matrix GreaterEqualTestScalar(Matrix a, long s)
{
  long nrows = rows(a);
  long ncols = cols(a);
  Matrix ret =  zeros(nrows, ncols);

  for (long i =0; i< nrows; i++)
    {
      for (long j =0; j< ncols; j++)
	{
	  if( (a[M(i, j, ncols)] >= (s-TOL)) )
	    {
	      ret[M(i, j, ncols)] = 1;
	    }
	}
    }
  return(ret);
} //end of GreaterEqualTestScalar

Matrix LessEqualTestScalar(Matrix a, double s)
{
  long nrows = rows(a);
  long ncols = cols(a);
  Matrix ret =  zeros(nrows, ncols);

  for (long i =0; i< nrows; i++)
    {
      for (long j =0; j< ncols; j++)
	{
	  if( (a[M(i, j, ncols)] <= (s+TOL)) )
	    {
	      ret[M(i, j, ncols)] = 1;
	    }
	}
    }
  return(ret);
} //end of LessEqualTestScalar

Matrix VectorAnd(Matrix a, Matrix b)
{
  long nrows = rows(a);
  Matrix ret =  zeros(nrows, 1);

  for (long i =0; i< nrows; i++)
    {
      if( (a[i] == 1) &  (b[i]== 1) )
	{
	  ret[i] = 1;
	}
    }
  return(ret);
} //end of VectorAnd

Matrix EqualityTestMatrix(Matrix a, Matrix s)
{
  long nrows = rows(a);
  long ncols = cols(a);
  Matrix ret =  zeros(nrows, ncols);

  long scols = cols(s);

  if (scols==1)
    {
      for (long i =0; i< nrows; i++)
	{
	  for (long j =0; j< ncols; j++)
	    {
	      if( (a[M(i, j, ncols)] < (s[i]+TOL)) &  (a[M(i, j, ncols)] > (s[i]-TOL)) )
		{
		  ret[M(i, j, ncols)] = 1;
		}
	    }
	}
    }
  else if (scols==ncols)
    {
      for (long i =0; i< nrows; i++)
	{
	  for (long j =0; j< ncols; j++)
	    {
	      if( (a[M(i, j, ncols)] < (s[M(i, j, ncols)]+TOL)) &  
		  (a[M(i, j, ncols)] > (s[M(i, j, ncols)]-TOL)) )
		{
		  ret[M(i, j, ncols)] = 1;
		}
	    }
	}
    }
  else 
    {
      printf("ASSERTION in EqualityTestMatrix\n");
    }

  return(ret);
} //end of EqualityTestMatrix

/* cumsum */
Matrix cumsum(Matrix a)
{
  long nrows = rows(a);
  Matrix ret = zeros(nrows, 1);

  ret[0] = a[0];
  for (long i = 1;  i < nrows; i++) 
    {
      ret[i] = ret[i-1] + a[i];
    }

  return ret;
} //end of cumsum

//! Calculate the sum of all of the elements in a  Matrix
double sum (const Matrix & A)
{
  double ret=0;
  long ncols = cols(A);
  long i;

  Matrix sumvec = sumc(A);

  for (i=0; i<ncols; i++)
    {
      ret = ret + sumvec[i];
    }

  return ret;
}

/*DISPLAY This function will display the double matrix stored in an mxArray.
 * This function assumes that the mxArray passed as input contains double
 * array.
 */
void display(Matrix A)
{
    int i=0, j=0; /* loop index variables */
    int r=0, c=0; /* variables to store the row and column length of the matrix */
    int count=0;

    /* Get the size of the matrix */
    r = rows(A);
    c = cols(A);
    
    /* Loop through the data and display the same in matrix format */
    for( i = 0; i < r; i++ ){
      for( j = 0; j < c; j++){
	printf("%4.2lf\t",A[count]);
	count++;
      }
      printf("\n");
    }
    printf("\n");
}
