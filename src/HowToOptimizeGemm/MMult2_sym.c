/* Create macros so that the matrices are stored in column-major order */

#define A(i,j) a[ (j)*lda + (i) ]
#define B(i,j) b[ (j)*ldb + (i) ]
#define C(i,j) c[ (j)*ldc + (i) ]

/* Routine for computing C = A * B + C */

void AddDot( int, double *, int, double *, double * );

/* m,n,k correspond to A[MxP], B[PxN], C[MxN] */
void MY_MMult( int m, int n, int p, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc )
{
  int i,k;  //really just awful.
  for(i=0; i < m/2; i++)  //top half of A
  {
    for(k=0; k < n; k++)  //every col of B
    {
      AddDot(p, &A(i,0), lda, &B(0,k), &C(i,k));  //this is causing a core dump for matrices over 200... :/
      C(m-i,n-k) = C(i,k);
    }
  }
  if(m%2 > 0)
  {
    for(i=0; i < n; i++)
    {
      AddDot(p, &A(m/2 + 1, 0), lda, &B(0,i), &C(m/2 + 1, i));
    }
  }
}
/* Create macro to let X( i ) equal the ith element of x */

#define X(i) x[ (i)*incx ]

void AddDot( int k, double *x, int incx,  double *y, double *gamma )
{
  /* compute gamma := x' * y + gamma with vectors x and y of length n.

     Here x starts at location x with increment (stride) incx and y starts at location y and has (implicit) stride of 1.
  */
 
  int p;

  for ( p=0; p<k; p++ ){
    *gamma += X( p ) * y[ p ];     
  }
}
