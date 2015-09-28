#include "heisen.h"

int main(int argc, char** argv) {
  HeisenCalculation calc(4);
  vector<double> e = calc.eigenvalues();
  print_vector(e);
  return 0;
}
