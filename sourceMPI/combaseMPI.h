/*
 *  combaseMPI.h
 *  Ulambator v0.3
 *
 *  Created by Mathias Nagel on 29.01.10.
 *  Copyright 2010 EPFL. All rights reserved.
 *
 */

/*
 *  LAPACK/BLAS header
 *  I use this header to declare LAPACK and BLAS routines.
 *	Furthermore I built a wrapper that changes from call by reference to call by value.
 *	However this is optional and up to the gusto of the user.
 *	Compile with -framework Accelerate or just -llapack -lblas
 *
 *
 *  Created by Mathias Nagel on 13.12.09.
 *  Copyright 2009 EPFL. All rights reserved.
 *
 */

int dgesv(double *A, double *b, int n);									// Matrix solver

int dsysv(double *A, double *b, int n);									// Matrix solver

int dgels(double *A, double *b, int m, int n, double *x);							// overdetermined Matrix solver

void dgemv(double *A, double *x, double *y, int m, int n, double alpha, double beta);									// Vektor multiplication

void dgemm(double *A, double *B, double *C, int m, int n, int k, double alpha, double beta);

void dg_ctof(double **in, double *lout, int rows, int cols);

int dgsev2rhs(double *A, double *b1, double *b2, int n);

int dgeinv(double *A, int N);

int schur(double *A, double *B, double *C, double *D, double *f1, double *f2, int M, int N);

int dgesch(double *A, double *f1, int M, int N);

extern "C" {
	extern void dgesv_(int *n, int *nrhs, double *a, int *lda,
					   int *ipiv, double *b, int *ldb, int *info);
    
    extern void dsysv_(char *uplo,int *n, int *nrhs, double *a, int *lda,
					   int *ipiv, double *b, int *ldb, double *work, int *lwork, int *info);
	
	extern void dgemv_(char *trans, int *m, int *n,
					   double *alpha, double *a, int *lda, double *x, int *incx,
					   double *beta,  double *y, int *incy ); 
	
	extern void dgels_(char *trans, int *m, int *n, int *nrhs,
					   double *a, int *lda, double *b, int *ldb, double *work, int *lwork,
					   int *info );
    extern void dgemm_(char *transA, char *transB, int*m, int*n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double* C, int *ldc);
    
    extern void dgetri_(int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info);
    extern void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv , int *info);
}


void bessk01(double x,double &k0,double &k1);											// Computation of Bessel functions
void bessi12(double x,double &i1,double &i2);											// Computation of Bessel functions

//Geometry Giacomo
void SplineSymmetric(double* x ,double* y, double *ax, double *bx, double *cx, double *dx, double *ay, double *by, double *cy, double *dy, int n);
double NormalVectorNxR(double& bx, double& by);
double NormalVectorNyR(double& bx, double& by);
double NormalVectorNxL(double& bx, double& cx, double& dx, double& by, double& cy, double& dy);
double NormalVectorNyL(double& bx, double& cx, double& dx, double& by, double& cy, double& dy);
double CurvSplineR(double &bx, double &cx, double &by, double &cy);
double CurvSplineL(double &bx, double &cx, double &dx, double &by, double &cy, double &dy);
void ChooseDistributionNodes(double*, double*, double*, double*, double, int, int, double, int, double);

//green functions Giacomo
void computeS(double,double,double,double,bool,bool,double&, double&, double&, double&);
void computeSSplines(double, double, double, double, double, double, double, double, bool, bool, bool, double &, double &, double &, double &);
void computeGTlinear1layer(double, double, double,double,double,double,double&, double&,double&,double&, double&, double&, double&, double&);
void computeGTlinearSPlines1layer(double, double, double,double, double, double, double, double, double, double, double, double, double&, double&,double&,double&, double&, double&,double&, double&);
double CenterMassAxis(double*,double*,int);
