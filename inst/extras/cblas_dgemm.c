/*
 *
 * cblas_dgemm.c
 * This program is a C interface to dgemm.
 * Written by Keita Teranishi
 * 4/8/1998
 *
 */

#include <R.h>
#include <R_ext/Applic.h> /* No longer contains R blas declarations */
#include <R_ext/BLAS.h> /* New R blas declarations */

#include "cblas.h"
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double  *A,
                 const int lda, const double  *B, const int ldb,
                 const double beta, double  *C, const int ldc)
{
   char TA, TB;   
#ifdef F77_CHAR
   F77_CHAR F77_TA, F77_TB;
#else
   #define F77_TA &TA  
   #define F77_TB &TB  
#endif

#ifdef F77_INT
   F77_INT F77_M=M, F77_N=N, F77_K=K, F77_lda=lda, F77_ldb=ldb;
   F77_INT F77_ldc=ldc;
#else
   #define F77_M M
   #define F77_N N
   #define F77_K K
   #define F77_lda lda
   #define F77_ldb ldb
   #define F77_ldc ldc
#endif

   /* the follow two vars were originally in Carla's.h as ex-terns, but that runs into conflicts on system with a
      preexisting cblas */
   int CBLAS_CallFromC;
   int RowMajorStrg;

   RowMajorStrg = 0;
   CBLAS_CallFromC = 1;

   if( Order == CblasColMajor )
   {
      if(TransA == CblasTrans) TA='T';
      else if ( TransA == CblasConjTrans ) TA='C';
      else if ( TransA == CblasNoTrans )   TA='N';
      else 
      {
	error("cblas_dgemm","Illegal TransA setting, %d\n", TransA);
	/* cblas_xerbla(2, "cblas_dgemm","Illegal TransA setting, %d\n", TransA); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      if(TransB == CblasTrans) TB='T';
      else if ( TransB == CblasConjTrans ) TB='C';
      else if ( TransB == CblasNoTrans )   TB='N';
      else 
      {
	error("cblas_dgemm","Illegal TransB setting, %d\n", TransB);
	/* cblas_xerbla(3, "cblas_dgemm","Illegal TransB setting, %d\n", TransB); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }

      #ifdef F77_CHAR
         F77_TA = C2F_CHAR(&TA);
         F77_TB = C2F_CHAR(&TB);
      #endif

	 F77_CALL(dgemm)(F77_TA, F77_TB, &F77_M, &F77_N, &F77_K, &alpha, A,
       &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc);
   } else if (Order == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if(TransA == CblasTrans) TB='T';
      else if ( TransA == CblasConjTrans ) TB='C';
      else if ( TransA == CblasNoTrans )   TB='N';
      else 
      {
	error("cblas_dgemm","Illegal TransA setting, %d\n", TransA);
	/* cblas_xerbla(2, "cblas_dgemm","Illegal TransA setting, %d\n", TransA); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      if(TransB == CblasTrans) TA='T';
      else if ( TransB == CblasConjTrans ) TA='C';
      else if ( TransB == CblasNoTrans )   TA='N';
      else 
      {
	error("cblas_dgemm","Illegal TransB setting, %d\n", TransB);
	/* cblas_xerbla(2, "cblas_dgemm","Illegal TransB setting, %d\n", TransB); */
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      #ifdef F77_CHAR
         F77_TA = C2F_CHAR(&TA);
         F77_TB = C2F_CHAR(&TB);
      #endif

	 F77_CALL(dgemm)(F77_TA, F77_TB, &F77_N, &F77_M, &F77_K, &alpha, B,
			 &F77_ldb, A, &F77_lda, &beta, C, &F77_ldc);
   } 
   else  {
     error("cblas_dgemm", "Illegal Order setting, %d\n", Order);
     /* cblas_xerbla(1, "cblas_dgemm", "Illegal Order setting, %d\n", Order); */
   }
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}


