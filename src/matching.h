#define M(ROW,COL,NCOLS) (((ROW)*(NCOLS))+(COL))
#define TOL 0.0000000001

/* Use CBLAS and Nate Optimizations */
#define __NBLAS__

// my function declarations
double sum (const Matrix & A);
double min_scalar (double a, double b);
double max_scalar (double a, double b);
Matrix multi_scalar (Matrix a, Matrix b);
Matrix col_assign (Matrix a, Matrix b, long col);
Matrix row_assign (Matrix a, Matrix b, long row);
Matrix diagCreate (Matrix a);
Matrix EqualityTestScalar(Matrix a, double s);
Matrix EqualityTestMatrix(Matrix a, Matrix s);
Matrix GreaterEqualTestScalar(Matrix a, long s);
Matrix LessEqualTestScalar(Matrix a, double s);
Matrix VectorAnd(Matrix a, Matrix b);
Matrix cumsum(Matrix a);
void display(Matrix A);





