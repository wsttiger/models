#include "heisen.h"

int main(int argc, char** argv) {
// //  for (int i = 0; i <= 4; i++) {
//     { int i = -1;
//     printf("%d\n", i);
//     HeisenCalculation calc(12,i);
//     vector<double> mat = calc.make_matrix();
//     int n = calc.get_nstates();
//     vector<double> e(n,0.0);
//     vector<double> ev(n*n, 0.0);
//     printf("\n");
//     print_matrix(mat,n,n);
//     diag_matrix(mat,n,e,ev);
//     print_vector(e);
// //    printf("Lanczos: lowest eigenvalue is %15.8f\n\n", e[0]);
//   }
  HeisenCalculation calc(12);
  vector<double> e = calc.eigenvalues();
  print_vector(e);
  return 0;
}
