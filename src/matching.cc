/*  
Jasjeet S. Sekhon <sekhon@berkeley.edu>
HTTP://sekhon.berkeley.edu/
UC Berkeley

2006/10/24
Under the GNU Public License Version 2
*/

/* We are using the include cblas header which will direct to the central R BLAS */

#include "scythematrix.h"

using namespace SCYTHE;
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

#if defined(__darwin__) || defined(__APPLE__)
#include <vecLib/cblas.h>
#else
#define INTERNAL_CBLAS
#include <R_ext/Applic.h>
#endif

#include "matching.h"

extern "C"
{
#include <Rdefines.h>

#ifdef INTERNAL_CBLAS
#include "cblas.h"

  void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
		   const int K, const double alpha, const double  *A,
		   const int lda, const double  *B, const int ldb,
		   const double beta, double  *C, const int ldc);
#endif

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                http://www.eso.org/~ndevilla/median/

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


  double kth_smallest(double *a, long n, long k)
  {
    long i,j,l,m ;
    double x, tmp;
    
    l=0 ; m=n-1 ;
    while (l<m) {
      x=a[k] ;
      i=l ;
      j=m ;
      do {
	while (a[i]<x) i++ ;
	while (x<a[j]) j-- ;
	if (i<=j) {
	  tmp=a[i];
	  a[i]=a[j];
	  a[j]=tmp;
	      i++ ; j-- ;
	}
      } while (i<=j) ;
      if (j<k) l=i ;
      if (k<i) m=j ;
    }
    return a[k] ;
  } // end of kth_smallest
  
  //Est Func
  SEXP EstFuncC (SEXP I_N, SEXP I_All, SEXP I_length, SEXP I_Y, SEXP I_Tr, SEXP I_weight, SEXP I_indx)
	      
  {
    SEXP ret;
    
    long N, All, i, l, j, k, SumFoo, length;
    double SumFooWeight, SumIndx3;
    long ic=3;

    N = asInteger(I_N);
    All = asInteger(I_All);
    length = asInteger(I_length);

    // In vars
    Matrix Tr = Matrix(N, 1);
    Matrix weight = Matrix(N, 1);
    Matrix indx = Matrix(length, ic);
    Matrix Y = Matrix(N, 1);
    
    // Extra vars
    Matrix Kcount = Matrix(N, 1);
    Matrix KKcount = Matrix(N, 1);
    Matrix YCAUS = Matrix(N, 1);
    Matrix foo = Matrix(N,1);
    Matrix WFW = Matrix(N,1);
    Matrix WWFW = Matrix(N, 1);

    k=0;
    //rows and colums are fliped!! j,i != i,j
    for(j=0;j<ic; j++)
      {
	for(i=0;i<length; i++)
	  {
	    indx[M(i,j,ic)] = REAL(I_indx)[k];
	    
	    // count from offest 0 instead of R's offset 1
	    if(j<2)
	      {
		indx[M(i,j,ic)] = indx[M(i,j,ic)] - 1;	
	      }
	    k++;
	  }
      }
    
    for(i=0;i<N; i++)
      {
	Y[i] = REAL(I_Y)[i];
	Tr[i] = REAL(I_Tr)[i];
	weight[i] = REAL(I_weight)[i];
      }
    
    
    //Master Loop
    for (i=0; i<N; i++)
      {
	if ( ((int) Tr[i]==1 & All!=1) | All==1)
	  {
	    
	    SumFooWeight = 0;
	    SumFoo = 0;
	    SumIndx3 = 0;
	    for (l=0; l<length; l++)
	      {
		if((int) indx[M(l,0,ic)]==i)
		  {
		    SumFoo++;
		    
		    k = (int) indx[M(l,1,ic)];
		    
		    foo[k] = 1.0;
		    
		    SumIndx3 = SumIndx3 + indx[M(l,2,ic)];
		    SumFooWeight = SumFooWeight + weight[k];
		  } // end of if
	      }//end l for
	  } //Tr[i]
	
	if (SumFoo > 0)
	  {
	    for (l=0; l<length; l++)
	      {
		
		if((int) indx[M(l,0,ic)]==i)
		  {
         	    k = (int) indx[M(l,1,ic)];
		    
		    WFW[k]  = weight[k]*foo[k]/SumFooWeight;
		    WWFW[k] = weight[k]*WFW[k]/SumFooWeight;
		    
		    Kcount[k] = Kcount[k] + weight[i] * WFW[k];
		    
		    KKcount[k] = KKcount[k] + weight[i]*WWFW[k];
		    
		    //if ( ((int) Tr[i]==1 & All!=1) | All==1)
		      YCAUS[i] = Y[k]*indx[M(l,2,ic)]/SumIndx3 + YCAUS[i];
		  } // endof if
	      } // l for loop
		
		if ((int) Tr[i]==1)
		  {
		    YCAUS[i] = Y[i] - YCAUS[i];
		  } 
		else 
		  {
		    YCAUS[i] = YCAUS[i] - Y[i];
		  }
	      } //if SumFoo
      }// end of master N loop

	    
    PROTECT(ret=allocMatrix(REALSXP, N, 3));
    
    // stack up cbind(YCAUS, Kcount, KKcount)
    k = 0;
    for( i = 0; i < N; i++ )
      {
	REAL(ret)[k] = YCAUS[i];
	k++;
      }
    for( i = 0; i < N; i++ )
      {
	REAL(ret)[k] = Kcount[i];
	k++;
      }
    for( i = 0; i < N; i++ )
      {
	REAL(ret)[k] = KKcount[i];
	k++;
      }
    UNPROTECT(1);
    return(ret); 
  } //end of EstFuncC

  SEXP FasterMatchC(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
		    SEXP I_ww, SEXP I_Tr, SEXP I_X)
  {
    SEXP ret;
    
    long N, xvars, All, M; 
    double cdd, Distmax, Wi_ratio;;
    
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
    
    Matrix INN = seqa(1, 1, N);
    // TT is just equal to INN

#ifdef __NBLAS__
    Matrix DX(N, xvars), ZX(N, xvars), Dist(N, xvars);
#else
    Matrix index_onesN = ones(N, 1);
    Matrix xx;
    Matrix DX, Dist;
#endif

    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0;
    
    Matrix DistPot, tt, xmat, I, IM, W;
    Matrix ACTMAT, POTMAT;

    // Temporary size for these matrices to avoid rbind
    int NM = N*M*100;
    Matrix tI  = Matrix(NM,1);
    Matrix tIM = Matrix(NM,1);
    Matrix tW  = Matrix(NM,1);
    
    //    tret = gettimeofday(&tv2,&tz2);
    //    long secs = tv2.tv_sec - tv1.tv_sec;
    //    long msecs = tv2.tv_usec - tv1.tv_usec;
    //    double actual = ((double) secs*1000000+ (double) msecs)/1000000;
    //    printf("actual: %lf, secs: %d, msecs: %d\n", actual, secs, msecs);
    
    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    //double *S = (double *) malloc(N*sizeof(double));  
    
    // struct timeval tv1, tv2;	// Required "timeval" structures (see man pg.)
    // struct timezone tz1, tz2;	// Required "timezone" structures (see man pg.)
    // long tret = gettimeofday(&tv1,&tz1);
    
    int MatchCount=0;
    int MCindx=0;
    int overFirstNM=1;
    for(i=0; i < N; i++)
    {
      // treatment indicator for observation to be matched        
      TREATi = (int) Tr[i];
      
      // proceed with all observations if All==1
      // but only with treated observations if All=0        
      if ( (TREATi==1 & All!=1) | All==1 )
      {
#ifdef __NBLAS__
        // this loop is equivalent to the matrix X less the matrix A, which is
        // the product of N rows by 1 column of 1.0 and row R of X.
        double *dest = ZX.data;
        double *src  = X.data;
        double *row  = X.data + (i * xvars);
        for (int jj = 0; jj < N; ++jj) {
          for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
            *dest = *src - row[kk];
          }
        }
        
        if (xvars>1)
	  {
	    //JSS
	    // do the second multiplication with dgemm, multiplying the matrix
	    // above, D, by the transpose of matrix W.

	    cblas_dgemm(CblasRowMajor,// column major
			CblasNoTrans, // A not transposed
			CblasTrans,   // B transposed
			xvars,        // M
			N,            // N
			xvars,        // K
			1.0,          // alpha, (alpha * A * B)
			ww.data,      // A
			xvars,        // leading dimension for A
			ZX.data,      // B
			xvars,        // leading dimension for B
			0.0,          // beta, (beta * C)
			DX.data,      // C
			N);           // leading dimension for C

          DX.multi_scalar(DX);
          std::swap(DX.colsize, DX.rowsize);
          
          Dist = sumc(DX);

          std::swap(Dist.colsize, Dist.rowsize); // transpose 1 x N -> N x 1
          std::swap(DX.colsize, DX.rowsize);
        } else {
          // do the second multiplication with dgemm, multiplying the matrix
          // above, D, by the transpose of matrix W.
          cblas_dgemm(CblasRowMajor, // column major
                      CblasNoTrans, // A not transposed
                      CblasTrans,   // B transposed
                      N,            // M
                      xvars,        // N
                      xvars,        // K
                      1.0,          // alpha, (alpha * A * B)
                      ZX.data,      // A
                      xvars,        // leading dimension for A
                      ww.data,      // B
                      xvars,        // leading dimension for B
                      0.0,          // beta, (beta * C)
                      DX.data,      // C
                      xvars);       // leading dimension for C
          
          DX.multi_scalar(DX);
          Dist = DX;
        } // end of xvars
#else
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
#endif /* end __NBLAS__ */
        
        // Dist distance to observation to be matched
        // is N by 1 vector	    
        
        // set of potential matches (all observations with other treatment)
        // JSS, note:logical vector
        POTMAT = EqualityTestScalar(Tr, 1-TREATi);
        
	// X's for potential matches
	DistPot = selif(Dist, POTMAT);
	//DistPot_size is a constant!! Fix it
	long DistPot_size = size(DistPot);
	
	//rsort_with_index(DistPot.data, order_DistPot, DistPot_size);
	//R_rsort(DistPot.data, DistPot_size);
	//rPsort(DistPot.data, DistPot_size, M);
	
	Distmax = kth_smallest(DistPot.data, DistPot_size, (M-1));        
        
        // selection of actual matches 
        // logical index
        ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
        ACTMAT = VectorAnd(POTMAT, ACTMAT);
        
        // distance to actual matches.  This is NEVER used.
        // ACTDIST = selif(Dist, ACTMAT);

	// counts how many times each observation is matched.
	double Mi = sum(ACTMAT);
	Wi_ratio = 1/Mi;

	// collect results
	MatchCount = MatchCount + (int) Mi;
	
	if(MatchCount > NM)
	  {
	    
	    NM = NM+N*M*100;
	    
	    if(overFirstNM > 0)
	      {
		printf("Increasing memory because of ties: allocating a matrix of size 3 times %d doubles.\n", NM);
		printf("I would be faster with the ties=FALSE option.\n");
		warning("Increasing memory because of ties.  I would be faster with the ties=FALSE option.");
	      }
	    
	    int OldMatchCount = MatchCount - (int) Mi;
	    
	    Matrix tI_tmp  = Matrix(NM,1);
	    Matrix tIM_tmp = Matrix(NM,1);
	    Matrix tW_tmp  = Matrix(NM,1);
	    
	    memcpy(tI_tmp.data, tI.data, OldMatchCount*sizeof(double));
	    memcpy(tIM_tmp.data, tIM.data, OldMatchCount*sizeof(double));
	    memcpy(tW_tmp.data, tW.data, OldMatchCount*sizeof(double));
	    
	    tI = tI_tmp;
	    tIM = tIM_tmp;
	    tW = tW_tmp;
	  }
	
	//foo1 = ones(ACTMATsum, 1)*(i+1);
	//memcpy(tI.data+MCindx, foo1.data, foo1.size*sizeof(double));
	//memcpy(tW.data+MCindx, Wi.data, Wi.size*sizeof(double));
	
	for (j=0; j < (int) Mi; j++)
	  {
	    tI.data[MCindx+j] = i+1;
	    tW.data[MCindx+j] = Wi_ratio;
	  }
	
	//foo1 = selif(INN, ACTMAT);
	//memcpy(tIM.data+MCindx, foo1.data, foo1.size*sizeof(double));
	//MCindx = MCindx+foo1.size;
	
	k=0;
	for (j=0; j<N; j++)
	  {
	    if(ACTMAT.data[j] > (1-TOL))
	      {
		tIM.data[MCindx+k] = j+1;
		k++;
	      }
	  }
	MCindx = MCindx+k;
      } // end of (TREATi==1 & All!=1) | All==1 )
    } //END OF i MASTER LOOP!

    // subset matrices to get back to the right dims
    if (MatchCount > 0)
      {
	I=Matrix(MatchCount, 1);
	IM=Matrix(MatchCount, 1);
	W=Matrix(MatchCount, 1);

	memcpy(I.data, tI.data, MatchCount*sizeof(double));
	memcpy(IM.data, tIM.data, MatchCount*sizeof(double));
	memcpy(W.data, tW.data, MatchCount*sizeof(double));
      }
    else
      {
	I=Matrix(1, 1);
	IM=Matrix(1, 1);
	W=Matrix(1, 1);
      }
    
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
    //free(S);
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
  } //end of FasterMatchC


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

#ifdef __NBLAS__
    Matrix DX(N, xvars), ZX(N, xvars);
#else
    Matrix index_onesN = ones(N, 1);
    Matrix xx;
    Matrix DX;
#endif

    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0, ACTMATsum = 0;
    
    Matrix Dist, DistPot, weightPot, tt, 
      weightPot_sort, weightPot_sum, Wi, xmat, I, IM, W;
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
#ifdef __NBLAS__
        // this loop is equivalent to the matrix X less the matrix A, which is
        // the product of N rows by 1 column of 1.0 and row R of X.
        double *dest = ZX.data;
        double *src  = X.data;
        double *row  = X.data + (i * xvars);
        for (int jj = 0; jj < N; ++jj) {
          for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
            *dest = *src - row[kk];
          }
        }
        
        if (xvars>1)
	  {
	    //JSS
	    // do the second multiplication with dgemm, multiplying the matrix
	    // above, D, by the transpose of matrix W.

	    cblas_dgemm(CblasRowMajor,// column major
			CblasNoTrans, // A not transposed
			CblasTrans,   // B transposed
			xvars,        // M
			N,            // N
			xvars,        // K
			1.0,          // alpha, (alpha * A * B)
			ww.data,      // A
			xvars,        // leading dimension for A
			ZX.data,      // B
			xvars,        // leading dimension for B
			0.0,          // beta, (beta * C)
			DX.data,      // C
			N);           // leading dimension for C

          DX.multi_scalar(DX);
          std::swap(DX.colsize, DX.rowsize);
          
          Dist = sumc(DX);

          std::swap(Dist.colsize, Dist.rowsize); // transpose 1 x N -> N x 1
          std::swap(DX.colsize, DX.rowsize);
        } else {
          // do the second multiplication with dgemm, multiplying the matrix
          // above, D, by the transpose of matrix W.
          cblas_dgemm(CblasRowMajor, // column major
                      CblasNoTrans, // A not transposed
                      CblasTrans,   // B transposed
                      N,            // M
                      xvars,        // N
                      xvars,        // K
                      1.0,          // alpha, (alpha * A * B)
                      ZX.data,      // A
                      xvars,        // leading dimension for A
                      ww.data,      // B
                      xvars,        // leading dimension for B
                      0.0,          // beta, (beta * C)
                      DX.data,      // C
                      xvars);       // leading dimension for C
          
          DX.multi_scalar(DX);
          Dist = DX;
        } // end of xvars
#else
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
#endif /* end __NBLAS__ */
        
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
        
        // distance to actual matches.  This is NEVER used.
        // ACTDIST = selif(Dist, ACTMAT);
        
        // counts how many times each observation is matched.
        double Mi = sum(multi_scalar(weight, ACTMAT));
        
#ifdef __NBLAS__
        foo1 = weight;
        foo1.multi_scalar(weight);
        foo1.multi_scalar(ACTMAT);
#else
        foo1 = multi_scalar(weight, weight);
        foo1 = multi_scalar(foo1, ACTMAT);
#endif
        
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
                  SEXP I_caliper, SEXP I_replace, SEXP I_ties,
                  SEXP I_ww, SEXP I_Tr, SEXP I_X, SEXP I_weight,
                  SEXP I_CaliperVec, SEXP I_Xorig,
                  SEXP I_restrict_trigger, SEXP I_restrict_nrow, SEXP I_restrict)
  {
    SEXP ret;
    
    long N, xvars, All, M, caliper, replace, ties, restrict_trigger, restrict_nrow, sum_caliper_drops=0,
      replace_count=0;
    double cdd, diff;
    
    long i, j, k, r, c;
    
    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    caliper = (long) asReal(I_caliper);
    replace = asInteger(I_replace);
    ties = asInteger(I_ties);
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

    // required for replace=0, to keep track of obs already used
    long *ReplaceVector;
    if(replace==0)
      {
	ReplaceVector = (long *) malloc(N*sizeof(long));
      }

    //Start random number generator if we are breaking ties
    if (ties==0)
      GetRNGstate();

    Matrix INN = seqa(1, 1, N);
    // TT is just equal to INN

#ifdef __NBLAS__
    Matrix DX(N, xvars), ZX(N, xvars);
#else
    Matrix index_onesN = ones(N, 1);
    Matrix xx;
    Matrix DX;
#endif

    Matrix Kcount  = zeros(N,1);
    Matrix KKcount = zeros(N,1);
    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0, ACTMATsum = 0;

    Matrix Dist, DistPot, weightPot, tt, 
      weightPot_sort, weightPot_sum, Wi, xmat, I, IM, IMt, W;
    Matrix ACTMAT, POTMAT(N, 1);

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
#ifdef __NBLAS__
	    // this loop is equivalent to the matrix X less the matrix A, which is
	    // the product of N rows by 1 column of 1.0 and row R of X.
	    double *dest = ZX.data;
	    double *src  = X.data;
	    double *row  = X.data + (i * xvars);
	    for (int jj = 0; jj < N; ++jj) {
	      for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		*dest = *src - row[kk];
	      }
	    }
	    
	    if (xvars>1)
	      {
		//JSS
		// do the second multiplication with dgemm, multiplying the matrix
		// above, D, by the transpose of matrix W.
		
		cblas_dgemm(CblasRowMajor,// column major
			    CblasNoTrans, // A not transposed
			    CblasTrans,   // B transposed
			    xvars,        // M
			    N,            // N
			    xvars,        // K
			    1.0,          // alpha, (alpha * A * B)
			    ww.data,      // A
			    xvars,        // leading dimension for A
			    ZX.data,      // B
			    xvars,        // leading dimension for B
			    0.0,          // beta, (beta * C)
			    DX.data,      // C
			    N);           // leading dimension for C
		
		DX.multi_scalar(DX);
		std::swap(DX.colsize, DX.rowsize);
		
		Dist = sumc(DX);
		
		std::swap(Dist.colsize, Dist.rowsize); // transpose 1 x N -> N x 1
		std::swap(DX.colsize, DX.rowsize);
	      } else {
	      // do the second multiplication with dgemm, multiplying the matrix
	      // above, D, by the transpose of matrix W.
	      cblas_dgemm(CblasRowMajor, // column major
			  CblasNoTrans, // A not transposed
			  CblasTrans,   // B transposed
			  N,            // M
			  xvars,        // N
			  xvars,        // K
			  1.0,          // alpha, (alpha * A * B)
			  ZX.data,      // A
			  xvars,        // leading dimension for A
			  ww.data,      // B
			  xvars,        // leading dimension for B
			  0.0,          // beta, (beta * C)
			  DX.data,      // C
			  xvars);       // leading dimension for C
	      
	      DX.multi_scalar(DX);
	      Dist = DX;
	    } // end of xvars
#else
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
#endif /* end __NBLAS__ */
	    
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

	    // Don't match observations we have already matched
	    if(replace_count > 0)
	      {
		for(j=0; j<replace_count; j++)
		  {
		    Dist[ ( ReplaceVector[j] - 1) ] = DOUBLE_XMAX;
		  }
	      } // end of replace_count

            // set of potential matches (all observations with other treatment)
            // JSS, note:logical vector
	    POTMAT = EqualityTestScalar(Tr, 1-TREATi);
	    
	    if (caliper==1)
	      {
		for (j=0; j<N; j++)
		  {
		    if((int) POTMAT[j]==1)
		      {
			for (k=0; k<xvars; k++)
			  {
			    diff = abs(Xorig[M(i, k, xvars)] - Xorig[M(j,k,xvars)]); 
			    if (diff > CaliperVec[k])
			      {
				Dist[j] = DOUBLE_XMAX;
				break;
			      }
			  }
		      }
		  } 
	      }//end of if caliper
	    
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

	    if (restrict_trigger==1 | caliper==1)
	      {
		if ( (Distmax+cdd) > DOUBLE_XMAX_CHECK)
		  {
		    sum_caliper_drops++;
		    continue;
		  }
	      } 

            // selection of actual matches 
            // logical index
	    ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	    ACTMAT = VectorAnd(POTMAT, ACTMAT);

	    if (ties==0)
	      {
		int Mii = (int) sum(ACTMAT);
		//Do we have ties?
		if (Mii > M)
		  {
		    IMt = selif(INN, ACTMAT);
		    int nties_broken = 0;
		    int ntiesToBreak = Mii - M;
		    while (nties_broken < ntiesToBreak) {
		      int idrop = (int) ( unif_rand()*(double) Mii);
		      k = (int) IMt[idrop];
		      if (k > (-1+TOL) )
			{
			  ACTMAT[k - 1] = 0;
			  IMt[idrop] = -1;
			  nties_broken++;
			} // end if
		    }// end of while loop
		  }
	      }// end of ties loop	    
	    
	    // distance to actual matches.  This is NEVER used.
	    // ACTDIST = selif(Dist, ACTMAT);
	    
            // counts how many times each observation is matched.
	    double Mi = sum(multi_scalar(weight, ACTMAT));
	    
            Kcount = Kcount + 
	      weight[i] * multi_scalar(weight, ACTMAT)/Mi;
	    
#ifdef __NBLAS__
	    foo1 = weight;
	    foo1.multi_scalar(weight);
	    foo1.multi_scalar(ACTMAT);
#else
	    foo1 = multi_scalar(weight, weight);
	    foo1 = multi_scalar(foo1, ACTMAT);
#endif
	    
	    KKcount = KKcount + (weight[i]* foo1)/(Mi*Mi);
	    
            Wi = selif(weight, ACTMAT);
	    Wi = weight[i]*Wi/Mi;
	    
	    ACTMATsum = (int) sumc(ACTMAT)[0];

	    //if no replacement
	    if(replace==0)
	      {
		IMt = selif(INN, ACTMAT);
		for (j=0; j<IMt.rowsize; j++)
		  {
		    ReplaceVector[replace_count] = (long) IMt.data[j];
		    replace_count++;
		  }
	      }//end of replace==0

	    // collect results
	    if (first==1)
	      {
		I = ones(ACTMATsum, 1)*(i+1);
		if(replace==0)
		  {
		    IM = IMt;
		  } 
		else 
		  {
		    IM = selif(INN, ACTMAT);
		  }
		W = Wi;
		first = 0;
	      }// end of first==1 
	    else 
	      {
		I = rbind(I, ones(ACTMATsum, 1)*(i+1));
		if(replace==0)
		  {
		    IM = rbind(IM, IMt);
		  } 
		else
		  {
		    IM = rbind(IM, selif(INN, ACTMAT));
		  }
		W = rbind(W, Wi);
	      } // end of i else
	    
	  } // end of (TREATi==1 & All!=1) | All==1 )
      } //END OF i MASTER LOOP!

    //Stop random number generator if we are breaking ties
    if (ties==0)
      PutRNGstate();
    
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
    if(replace==0)
      {
	free(ReplaceVector);
      }
    
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

  SEXP MatchLoopCfast(SEXP I_N, SEXP I_xvars, SEXP I_All, SEXP I_M, SEXP I_cdd,
		      SEXP I_caliper, SEXP I_replace, SEXP I_ties,
		      SEXP I_ww, SEXP I_Tr, SEXP I_X, 
		      SEXP I_CaliperVec, SEXP I_Xorig,
		      SEXP I_restrict_trigger, SEXP I_restrict_nrow, SEXP I_restrict)
  {
    SEXP ret;
    
    long N, xvars, All, M, caliper, replace, ties, restrict_trigger, restrict_nrow, sum_caliper_drops=0,
      replace_count=0;
    double cdd, diff, Distmax, Wi_ratio;
    
    long i, j, k, r, c;
    
    N = asInteger(I_N);
    xvars = asInteger(I_xvars);
    All = asInteger(I_All);
    M = asInteger(I_M);
    cdd = asReal(I_cdd);
    caliper = (long) asReal(I_caliper);
    replace = asInteger(I_replace);
    ties = asInteger(I_ties);
    restrict_nrow = asInteger(I_restrict_nrow);
    restrict_trigger = asInteger(I_restrict_trigger);
    
    //    struct timeval tv1, tv2;	// Required "timeval" structures (see man pg.)
    //    struct timezone tz1, tz2;	// Required "timezone" structures (see man pg.)
    //    long tret = gettimeofday(&tv1,&tz1);
    
    Matrix ww = Matrix(xvars, xvars);
    Matrix Tr = Matrix(N, 1);
    Matrix X = Matrix(N, xvars);
    
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

    // required for replace=0, to keep track of obs already used
    long *ReplaceVector;
    if(replace==0)
      {
	ReplaceVector = (long *) malloc(N*sizeof(long));
      }

    //Start random number generator if we are breaking ties
    if (ties==0)
      GetRNGstate();

    Matrix INN = seqa(1, 1, N);
    // TT is just equal to INN

#ifdef __NBLAS__
    Matrix DX(N, xvars), ZX(N, xvars), Dist(N, xvars);
#else
    Matrix index_onesN = ones(N, 1);
    Matrix xx;
    Matrix DX, Dist;
#endif

    Matrix foo1;
    Matrix foo2;

    /* set misc */
    int TREATi = 0;

    Matrix tt, xmat, I, IM,  IMt, W;
    Matrix ACTMAT, POTMAT(N, 1);

    // Temporary size for these matrices to avoid rbind
    int NM = N*M*100;
    if (ties==0)
      NM = N*M;

    Matrix tI  = Matrix(NM,1);
    Matrix tIM = Matrix(NM,1);
    Matrix tW  = Matrix(NM,1);

    //    tret = gettimeofday(&tv2,&tz2);
    //    long secs = tv2.tv_sec - tv1.tv_sec;
    //    long msecs = tv2.tv_usec - tv1.tv_usec;
    //    double actual = ((double) secs*1000000+ (double) msecs)/1000000;
    //    printf("actual: %lf, secs: %d, msecs: %d\n", actual, secs, msecs);

    //These are larger than needed; it is just easier this way
    int *order_DistPot = (int *) malloc(N*sizeof(int));  
    //double *S = (double *) malloc(N*sizeof(double));  
    Matrix DistPot=Matrix(N, 1);

    // struct timeval tv1, tv2;	// Required "timeval" structures (see man pg.)
    // struct timezone tz1, tz2;	// Required "timezone" structures (see man pg.)
    // long tret = gettimeofday(&tv1,&tz1);
    
    int MatchCount=0;
    int MCindx=0;
    int overFirstNM=1;
    for(i=0; i < N; i++)
      {
	// treatment indicator for observation to be matched        
	TREATi = (int) Tr[i];
	
	// proceed with all observations if All==1
	// but only with treated observations if All=0        
	if ( (TREATi==1 & All!=1) | All==1 )
	  {
#ifdef __NBLAS__
	    // this loop is equivalent to the matrix X less the matrix A, which is
	    // the product of N rows by 1 column of 1.0 and row R of X.
	    double *dest = ZX.data;
	    double *src  = X.data;
	    double *row  = X.data + (i * xvars);
	    for (int jj = 0; jj < N; ++jj) {
	      for (int kk = 0; kk < xvars; ++kk, ++dest, ++src) {
		*dest = *src - row[kk];
	      }
	    }

	    if (xvars>1)
	      {
		//JSS
		// do the second multiplication with dgemm, multiplying the matrix
		// above, D, by the transpose of matrix W.
		
		cblas_dgemm(CblasRowMajor,// column major
			    CblasNoTrans, // A not transposed
			    CblasTrans,   // B transposed
			    xvars,        // M
			    N,            // N
			    xvars,        // K
			    1.0,          // alpha, (alpha * A * B)
			    ww.data,      // A
			    xvars,        // leading dimension for A
			    ZX.data,      // B
			    xvars,        // leading dimension for B
			    0.0,          // beta, (beta * C)
			    DX.data,      // C
			    N);           // leading dimension for C
		
		DX.multi_scalar(DX);
		std::swap(DX.colsize, DX.rowsize);
		
		Dist = sumc(DX);
		
		std::swap(Dist.colsize, Dist.rowsize); // transpose 1 x N -> N x 1
		std::swap(DX.colsize, DX.rowsize);
	      } else {
	      // do the second multiplication with dgemm, multiplying the matrix
	      // above, D, by the transpose of matrix W.
	      cblas_dgemm(CblasRowMajor, // column major
			  CblasNoTrans, // A not transposed
			  CblasTrans,   // B transposed
			  N,            // M
			  xvars,        // N
			  xvars,        // K
			  1.0,          // alpha, (alpha * A * B)
			  ZX.data,      // A
			  xvars,        // leading dimension for A
			  ww.data,      // B
			  xvars,        // leading dimension for B
			  0.0,          // beta, (beta * C)
			  DX.data,      // C
			  xvars);       // leading dimension for C

	      DX.multi_scalar(DX);
	      Dist = DX;
	    } // end of xvars
#else
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
#endif /* end __NBLAS__ */
	    
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

	    // Don't match observations we have already matched
	    if(replace_count > 0)
	      {
		for(j=0; j<replace_count; j++)
		  {
		    Dist[ ( ReplaceVector[j] - 1) ] = DOUBLE_XMAX;
		  }
	      } // end of replace_count

            // set of potential matches (all observations with other treatment)
            // JSS, note:logical vector
	    POTMAT = EqualityTestScalar(Tr, 1-TREATi);

	    if (caliper==1)
	      {
		for (j=0; j<N; j++)
		  {
		    if((int) POTMAT[j]==1)
		      {
			for (k=0; k<xvars; k++)
			  {
			    diff = abs(Xorig[M(i, k, xvars)] - Xorig[M(j,k,xvars)]); 
			    if (diff > CaliperVec[k])
			      {
				Dist[j] = DOUBLE_XMAX;
				break;
			      }
			  }
		      }
		  } 
	      }//end of if caliper

            // X's for potential matches
            DistPot = selif(Dist, POTMAT);
	    //DistPot_size is a constant!! Fix it
	    long DistPot_size = size(DistPot);

	    //rsort_with_index(DistPot.data, order_DistPot, DistPot_size);
	    //R_rsort(DistPot.data, DistPot_size);
	    //rPsort(DistPot.data, DistPot_size, M);

	    Distmax = kth_smallest(DistPot.data, DistPot_size, (M-1));

	    if (restrict_trigger==1 | caliper==1)
	      {
		if ( (Distmax+cdd) > DOUBLE_XMAX_CHECK)
		  {
		    sum_caliper_drops++;
		    continue;
		  }
	      } 

            // selection of actual matches 
            // logical index
	    ACTMAT = LessEqualTestScalar(Dist,  (Distmax+cdd));
	    ACTMAT = VectorAnd(POTMAT, ACTMAT);

	    if (ties==0)
	      {
		int Mii = (int) sum(ACTMAT);
		//Do we have ties?
		if (Mii > M)
		  {
		    IMt = selif(INN, ACTMAT);
		    int nties_broken = 0;
		    int ntiesToBreak = Mii - M;
		    while (nties_broken < ntiesToBreak) {
		      int idrop = (int) ( unif_rand()*(double) Mii);
		      k = (int) IMt[idrop];
		      if (k > (-1+TOL) )
			{
			  ACTMAT[k - 1] = 0;
			  IMt[idrop] = -1;
			  nties_broken++;
			} // end if
		    }// end of while loop
		  }
	      }// end of ties loop	    
	    
	    // distance to actual matches.  This is NEVER used.
	    // ACTDIST = selif(Dist, ACTMAT);
	    
            // counts how many times each observation is matched.
	    double Mi = sum(ACTMAT);
	    //ACTMATsum = (int) sumc(ACTMAT)[0];

	    Wi_ratio = 1/Mi;
	    //Wi = ones(ACTMATsum, 1)*1/Mi;

	    //if no replacement
	    if(replace==0)
	      {
		IMt = selif(INN, ACTMAT);
		for (j=0; j<IMt.rowsize; j++)
		  {
		    ReplaceVector[replace_count] = (long) IMt.data[j];
		    replace_count++;
		  }
	      }//end of replace==0

	    // collect results
	    MatchCount = MatchCount + (int) Mi;

	    /*
	    if (MatchCount > NM)
	      {
		MatchCount = MatchCount - (int) Mi;
		continue;
	      }
	    */

	    if(MatchCount > NM)
	      {
		
		NM = NM+N*M*100;

		if(overFirstNM > 0)
		  {
		    printf("Increasing memory because of ties: allocating a matrix of size 3 times %d doubles.\n", NM);
		    printf("I would be faster with the ties=FALSE option.\n");
		    warning("Increasing memory because of ties.  I would be faster with the ties=FALSE option.");
		  }
		else 
		  {
		    printf("Increasing memory because of ties: allocating a matrix of size 3 times %d doubles.", NM);
		  }
		
		int OldMatchCount = MatchCount - (int) Mi;
		
		Matrix tI_tmp  = Matrix(NM,1);
		Matrix tIM_tmp = Matrix(NM,1);
		Matrix tW_tmp  = Matrix(NM,1);
		
		memcpy(tI_tmp.data, tI.data, OldMatchCount*sizeof(double));
		memcpy(tIM_tmp.data, tIM.data, OldMatchCount*sizeof(double));
		memcpy(tW_tmp.data, tW.data, OldMatchCount*sizeof(double));
		
		tI = tI_tmp;
		tIM = tIM_tmp;
		tW = tW_tmp;

		/* free(tI_tmp.data);
		free(tIM_tmp.data);
		free(tW_tmp.data); */
	      }

	    //foo1 = ones(ACTMATsum, 1)*(i+1);
	    //memcpy(tI.data+MCindx, foo1.data, foo1.size*sizeof(double));
	    //memcpy(tW.data+MCindx, Wi.data, Wi.size*sizeof(double));

	    for (j=0; j < (int) Mi; j++)
	      {
		tI.data[MCindx+j] = i+1;
		tW.data[MCindx+j] = Wi_ratio;
	      }
	    
	    if(replace==0)
	      {
		memcpy(tIM.data+MCindx, IMt.data, IMt.size*sizeof(double));
		MCindx = MCindx+IMt.size;
	      }
	    else
	      {
		//foo1 = selif(INN, ACTMAT);
		//memcpy(tIM.data+MCindx, foo1.data, foo1.size*sizeof(double));
		//MCindx = MCindx+foo1.size;
		
		k=0;
		for (j=0; j<N; j++)
		  {
		    if(ACTMAT.data[j] > (1-TOL))
		      {
			tIM.data[MCindx+k] = j+1;
			k++;
		      }
		  }
		MCindx = MCindx+k;
	      }
	  } // end of (TREATi==1 & All!=1) | All==1 )
      } //END OF i MASTER LOOP!

    //Stop random number generator if we are breaking ties
    if (ties==0)
      PutRNGstate();

    // subset matrices to get back to the right dims
    if (MatchCount > 0)
      {
	I=Matrix(MatchCount, 1);
	IM=Matrix(MatchCount, 1);
	W=Matrix(MatchCount, 1);

	memcpy(I.data, tI.data, MatchCount*sizeof(double));
	memcpy(IM.data, tIM.data, MatchCount*sizeof(double));
	memcpy(W.data, tW.data, MatchCount*sizeof(double));
      }
    else
      {
	I=Matrix(1, 1);
	IM=Matrix(1, 1);
	W=Matrix(1, 1);
      }

#if defined(NERVERDEF)
    printf("I.rowsize %d\n", I.rowsize);
    printf("tI.rowsize %d\n", tI.rowsize);

    printf("IM.rowsize %d\n", IM.rowsize);
    printf("tIM.rowsize %d\n", tIM.rowsize);

    printf("W.rowsize %d\n", W.rowsize);
    printf("tW.rowsize %d\n", tW.rowsize);

    printf("NM: %d\n", NM);
    printf("display\n");
    fflush(stdout);

    for(i=0; i<I.rowsize; i++)
      {
	printf("%d %d: %lf %lf %lf\n", i, i+1, 
	       I[i], 
	       IM[i],
	       W[i]); 
    	fflush(stdout);
      }
#endif

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
    //free(S);
    free(order_DistPot);
    if(replace==0)
      {
	free(ReplaceVector);
      }
    
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
  } //end of MatchLoopCfast

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

#ifdef __NBLAS__
Matrix multi_scalar (Matrix a, Matrix b)
{
  a.multi_scalar(b);
  return a;
} // multi_scalar
#else
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
#endif

#ifdef __NBLAS__
Matrix EqualityTestScalar(Matrix a, double s)
{
  for (long i = 0; i < a.size; ++i)
    a.data[i] = (a.data[i] < (s+TOL)) && (a.data[i] > (s-TOL)) ? 1 : 0;
  return a;
} //end of EqualityTestScalar
#else
Matrix EqualityTestScalar(Matrix a, double s)
{
  long nrows = rows(a);
  long ncols = cols(a);
  Matrix ret =  zeros(nrows, ncols);

  for (long i =0; i< nrows; i++)
    {
      for (long j =0; j< ncols; j++)
	{
	  if( (a[M(i, j, ncols)] < (s+TOL)) && (a[M(i, j, ncols)] > (s-TOL)) )
	    {
	      ret[M(i, j, ncols)] = 1;
	    }
	}
    }
  return(ret);
} //end of EqualityTestScalar
#endif

#ifdef __NBLAS__
Matrix GreaterEqualTestScalar(Matrix a, long s)
{
  for (long i = 0; i < a.size; ++i)
    a.data[i] = (a.data[i] >= (s-TOL)) ? 1 : 0;
  return a;
} //end of GreaterEqualTestScalar
#else
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
#endif

#ifdef __NBLAS__
Matrix LessEqualTestScalar(Matrix a, double s)
{
  for (long i = 0; i < a.size; ++i)
    a.data[i] = (a.data[i] <= (s+TOL)) ? 1 : 0;
  return a;
} //end of LessEqualTestScalar
#else
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
#endif

Matrix VectorAnd(Matrix a, Matrix b)
{
  long nrows = rows(a);
  Matrix ret =  zeros(nrows, 1);

  for (long i =0; i< nrows; i++)
    {
      if( (a[i] == 1) &&  (b[i]== 1) )
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


