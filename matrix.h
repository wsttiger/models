/*
 * matrix.h
 *
 *  Created on: Aug 22, 2011
 *      Author: wsttiger
 */

#include <vector>
#include <complex>

using std::vector;
using std::complex;

#ifndef MATRIX_H_
#define MATRIX_H_

extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a,
           int *lda, double *w, double *work, int *lwork,
           int *info);
extern "C" void dgetri(int* n, double* a, int *lda, int *ipiv,
            double* work, int* lwork, int* info);
extern "C" void dgetrf(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
extern "C" void zgetri(int* n, complex<double>* a, int *lda, int *ipiv,
            complex<double>* work, int* lwork, int* info);
extern "C" void zgetrf(int* m, int* n, complex<double>* a, int* lda, int* ipiv, int* info);

//*****************************************************************************
long fact(int n)
{
  long prod = 1;
  for (long i = 1; i <= n; i++)
  {
    prod *= i;
  }
  return prod;
}
//*****************************************************************************

//*****************************************************************************
int nchoosek(int n, int k)
{
  long rval = 1;
  for (long i = k+1; i <= n; i++)
  {
    rval *= i;
  }
  return (rval/fact(n-k));
}
//*****************************************************************************

//*****************************************************************************
double dotvec(const vector<double>& v1, const vector<double>& v2)
{
  double rval = 0.0;
  for (unsigned int iv = 0; iv < v1.size(); iv++)
  {
    rval += v1[iv]*v2[iv];
  }
  return rval;
}
//*****************************************************************************

//r = a*x + b*y
//*****************************************************************************
vector<double> wst_daxpy(const double& a, const vector<double>& x,
                          const double& b, const vector<double>& y)
{
  vector<double> rv(x.size(),0.0);
  for (unsigned int i = 0; i < rv.size(); i++)
  {
    rv[i] = a*x[i]+b*y[i];
  }
  return rv;
}
//*****************************************************************************

// in-place version of daxpy
// x = a*x + b*y
//*****************************************************************************
void wst_daxpy_inplace(const double& a, vector<double>& x,
                          const double& b, const vector<double>& y)
{
  for (unsigned int i = 0; i < x.size(); i++)
  {
    x[i] = a*x[i]+b*y[i];
  }
}
//*****************************************************************************

//*****************************************************************************
double norm2(const vector<double>& v)
{
  double rnorm = 0.0;
  for (unsigned int iv = 0; iv < v.size(); iv++)
  {
    rnorm += v[iv]*v[iv];
  }
  return sqrt(rnorm);
}
//*****************************************************************************

//*****************************************************************************
vector<double> matrix_vector_mult(vector<double> mat, vector<double> v)
{
  unsigned int szmat = mat.size();
  unsigned int sz2 = v.size();
  unsigned int sz1 = szmat / sz1;
  vector<double> rv(sz2,0.0);
  for (unsigned int i = 0; i < sz2; i++)
  {
    for (unsigned int j = 0; j < sz1; j++)
    {
      rv[i] += mat[i*sz1+j];
    }
  }
  return rv;
}
//*****************************************************************************

//*****************************************************************************
void normalize(vector<double>& v)
{
  double rnorm = 0.0;
  for (unsigned int iv = 0; iv < v.size(); iv++)
  {
    rnorm += v[iv]*v[iv];
  }
  rnorm = sqrt(rnorm);
  for (unsigned int iv = 0; iv < v.size(); iv++)
  {
    v[iv] /= rnorm;
  }
}
//*****************************************************************************

//*****************************************************************************
void print_vector(const vector<double>& v)
{
  for (unsigned int i = 0; i < v.size(); i++)
  {
    printf("  %10.5f\n",v[i]);
  }
}
//*****************************************************************************

//*****************************************************************************
void print_matrix(const vector<double>& mat, int m, int n)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%8.2f",mat[i*n+j]);
    }
    printf("\n");
  }
}
//*****************************************************************************

//*****************************************************************************
void print_matrix(const vector< complex<double> >& mat, int m, int n)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("   (%8.4f,%8.4f)",mat[i*n+j].real(),mat[i*n+j].imag());
    }
    printf("\n");
  }
}
//*****************************************************************************

//void diag_matrix2(const vector<double>& mat, int n,
//                  vector<double>& e, vector <double>& evec)
//{
//  MatrixXcf A = MatrixXcf::Zero(n,n);
//  for (int i = 0; i < n; i++)
//  {
//    for (int j = 0; j < n; j++)
//    {
//      A(i,j) = complex<double>(mat[i*n+j],0.0);
//    }
//  }
//
//  ComplexEigenSolver<MatrixXcf> ces;
//  ces.compute(A);
//  VectorXcf et = ces.eigenvalues();
//  MatrixXcf evt = ces.eigenvectors();
//
//  for (int i = 0; i < n; i++)
//  {
//    e[i] = real(et(i));
//    for (int j = 0; j < n; j++)
//    {
//      evec[i*n+j] = real(evt(i,j));
//    }
//  }
//
//}

//*****************************************************************************
vector<double> get_col_from_matrix(const vector<double> mat, int col, int m, int n)
{
  vector<double> rv(m,0.0);
  for (int i = 0; i < m; i++)
  {
    rv[i] = mat[i*n+col];
  }
  return rv;
}
//*****************************************************************************

//*****************************************************************************
void diag_matrix(const vector<double>& mat, int n,
                  vector<double>& e, vector <double>& evec)
{
  char jobz = 'V';
  char uplo = 'U';
  int info;
  double* a = new double[n*n];
  int lda = n;
  int lwork = 3*n-1;
  double *work = new double[lwork];
  double *et = new double[n];
  for (int i = 0; i < n*n; i++) a[i] = mat[i];

  dsyev_(&jobz, &uplo, &n, a, &lda, et, work, &lwork, &info);

  if (info != 0)
  {
    printf("[[Error:]] lapack::dsyev failed --- info = %d\n\n", info);
  }
  else
  {
    //printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++)
    {
      e[i] = et[i];
    }
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        evec[j*n+i] = a[i*n+j];
      }
    }

}
  delete a;
  delete work;
  delete et;
}
//*****************************************************************************

//*****************************************************************************
vector<double> invert_matrix(const vector<double>& mat, int n)
{
  vector<double> rv(n*n,0.0);
  int info;
  double* a = new double[n*n];
  int lda = n;
  int lwork = 3*n-1;
  int* ipiv = new int[n];
  for (int i = 0; i < n; i++) ipiv[i] = 0;
  double *work = new double[lwork];
  for (int i = 0; i < n*n; i++) a[i] = mat[i];

  dgetrf(&n, &n, a, &lda, ipiv, &info);
  if (info != 0)
  {
    printf("[[Error:]] lapack::dgetrf failed --- info = %d\n\n", info);
  }
  dgetri(&n, a, &lda, ipiv, work, &lwork, &info);

  if (info != 0)
  {
    printf("[[Error:]] lapack::dgetri failed --- info = %d\n\n", info);
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        rv[i*n+j] = a[i*n+j];
      }
    }

  }
  delete a;
  delete work;
  delete ipiv;

  return rv;
}
//*****************************************************************************

//*****************************************************************************
vector< complex<double> > invert_matrix(const vector< complex<double> >& mat, int n)
{
  vector< complex<double> > rv(n*n,0.0);
  int info;
  complex<double>* a = new complex<double>[n*n];
  int lda = n;
  int lwork = 3*n-1;
  int* ipiv = new int[n];
  for (int i = 0; i < n; i++) ipiv[i] = 0;
  complex<double>* work = new complex<double>[lwork];
  for (int i = 0; i < n*n; i++) a[i] = mat[i];

  zgetrf(&n, &n, a, &lda, ipiv, &info);
  if (info != 0)
  {
    printf("[[Error:]] lapack::zgetrf failed --- info = %d\n\n", info);
  }
  zgetri(&n, a, &lda, ipiv, work, &lwork, &info);

  if (info != 0)
  {
    printf("[[Error:]] lapack::zgetri failed --- info = %d\n\n", info);
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        rv[i*n+j] = a[i*n+j];
      }
    }

  }
  delete a;
  delete work;
  delete ipiv;

  return rv;
}
//*****************************************************************************

//*****************************************************************************
vector<double> random_vector(const int& nsize)
{
  vector<double> v(nsize,0.0);
  for (int i = 0; i < nsize; i++)
  {
    int i1 = rand();
    double t1 = (i1 % 100000)/100000.0;
    v[i] = t1;
  }
  normalize(v);
  return v;
}
//*****************************************************************************

//*****************************************************************************
bool is_equals(const vector<double>& a, const vector<double>& b, double tol = 1e-10)
{
  for (unsigned int i = 0; i < a.size(); i++)
  {
    if (fabs(a[i]-b[i]) > tol) return false;
  }
  return true;
}
//*****************************************************************************

#endif /* MATRIX_H_ */
