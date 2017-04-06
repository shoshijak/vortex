
#include "simulation.h"
#include <iostream>

int main()
{
  std::cout << "Starting VORTEX application - with OpenMP\n";

  double theta = 0.5;
  bool verify = true;
  double tol_ = 1e-8; // is supposed to be the same as tol in run_tests.cpp, but I was having linker errors

  double extT, mrtT, srtT, reoT, bldT, evaT, potT;
  int n, nnodes, NDST;
  int numtests = 5;
  bool printeach = true;

  run_test(extT, mrtT, srtT, reoT, bldT, evaT, potT,
           n, nnodes, NDST, theta, tol_, verify,
           numtests, printeach);

  return 0;
}

