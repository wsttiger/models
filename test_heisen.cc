#include "heisen.h"

int main(int argc, char** argv) {
  ifstream f("read.in");
  int nsites = 4;
  int sector = -1;
  if (f.is_open()) {
    f >> nsites; 
    f >> sector; 
    printf("read.in --> nsites = %d     sector = %d\n", nsites, sector);
    f.close();
  }
  HeisenCalculation calc(nsites,sector);
//  const auto tstart = std::chrono::system_clock::now();
  vector<double> e = calc.eigenvalues();
//  const auto tstop = std::chrono::system_clock::now();
//  const std::chrono::duration<double> time_elapsed = tstop - tstart;
  print_vector(e);
//  std::cout << "Calculation took " << time_elapsed.count() << " s" << std::endl;
  return 0;
}
