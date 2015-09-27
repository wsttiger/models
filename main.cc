/*
 * main.cc
 *
 *  Created on: Apr 20, 2011
 *      Author: wsttiger
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <bitset>
#include <vector>
#include <map>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "matrix.h"
#include "hubbard.h"

using std::cout;
using std::endl;
using std::bitset;
using std::vector;
using std::map;
using std::pair;
using std::string;
using std::complex;

using namespace Eigen;

#define NSIZE 26



//void test_eigen3()
//{
//  MatrixXcf A = MatrixXcf::Random(4,4);
//  cout << "Here is a random 4x4 matrix, A:" << endl << A << endl << endl;
//
//  ComplexEigenSolver<MatrixXcf> ces;
//  ces.compute(A);
//  VectorXcf ev = ces.eigenvalues();
//  MatrixXcf evv = ces.eigenvectors();
//
//  printf("WST: the lowest eigenvalue: %15.7f\n\n", real(ev[0]));
//
//  cout << "The eigenvalues of A are:" << endl << ces.eigenvalues() << endl;
//  cout << "The matrix of eigenvectors, V, is:" << endl << ces.eigenvectors() << endl << endl;
//
//  complex<float> lambda = ces.eigenvalues()[0];
//  cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
//  VectorXcf v = ces.eigenvectors().col(0);
//  cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
//  cout << "... and A * v = " << endl << A * v << endl << endl;
//
//  cout << "Finally, V * D * V^(-1) = " << endl
//       << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << endl;
//}

//void test_diag()
//{
//  int n = 3;
//  vector<double> mat(n,0.0);
//
//  for (int i = 0; i < n; i++)
//  {
//    for (int j = 0; j < n; j++)
//    {
//      mat[i*n+j] = i*n+j;
//    }
//  }
//
//  print_matrix(mat,n,n);
//  vector<double> e(n,0.0);
//  vector<double> ev(n*n,0.0);
//}

void hubbard_w_matrix()
{
  // create HubbardCalculation object
  HubbardCalculation hc(3,2,3,3,8.0,1.0);
  printf("number of states: %d\n\n", hc.get_nstates());

  int nstates = hc.get_nstates();
  vector<double> hmat = hc.make_matrix(hc.get_states());
  vector<double> e(nstates,0.0);
  vector<double> ev(nstates*nstates,0.0);
  diag_matrix(hmat,nstates,e,ev);
  printf("diag matrix: %15.8f\n\n", e[0]);
}

void hubbard_w_lanczos()
{
  // create HubbardCalculation object
  HubbardCalculation hc(3,2,3,3,8.0,1.0);
  printf("number of states: %d\n\n", hc.get_nstates());

  Lanczos<HubbardCalculation> l(&hc,100);
  l.run();

  vector<double> ev1 = l.lowstate();
  //print_vector(ev1);
  hc.compute_1p_greens_function(-10.0,10.0,200,0.1);
}

void hubbard_w_lanczos_from_file()
{
  // create HubbardCalculation object
  HubbardCalculation hc("hubbard.in");
  printf("number of states: %d\n\n", hc.get_nstates());

  //vector<double> mat = hc.make_matrix();
  //print_matrix(mat,hc.get_nstates(), hc.get_nstates());
  //exit(EXIT_FAILURE);

  Lanczos<HubbardCalculation> l(&hc,100);
  l.run();

  hc.make_matrix();

//  vector<double> ev1 = l.lowstate();
//  //print_vector(ev1);
////  hc.compute_1p_greens_function(-32.0,16.0,800,0.1);
//  hc.compute_1p_greens_function(-3.0,6.0,8,0.1);
//  hc.compute_1p_greens_function_matrix(-3.0,6.0,8,0.1);
}

void test_hubbard_w_openmp()
{
  // create HubbardCalculation object
  HubbardCalculation hc(4,4,6,6,8.0,1.0);
  unsigned int nstates = hc.get_nstates();
  printf("number of states: %d\n\n", nstates);

  // create random vector
  vector<double> v1 = random_vector(nstates);
  // apply hubbard hamiltonian to v1 without openmp
  vector<double> rv1 = hc.apply(v1);

  unsigned int ntimes = 200;
  for (unsigned int i = 0; i < ntimes; i++)
  {
    vector<double> rv2 = hc.apply_w_openmp(v1);
    string result = (is_equals(rv1,rv2)) ? "PASS!" : "FAIL!";
    cout << "trial: " << i << "     " << result << endl;
  }
}

int main(int argc, char** argv)
{
  hubbard_w_lanczos_from_file();

//  int tid;
//  printf("Hello world from threads:\n");
//  #pragma omp parallel private (tid)
//  {
//    tid = omp_get_thread_num();
//    printf("Hello from thread: <%d>\n", tid);
//  }
//  printf("I am sequential now.\n");
//
//  test_hubbard_w_openmp();
//  return 0;

  return 0;
}

