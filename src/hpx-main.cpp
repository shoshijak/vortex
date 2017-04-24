
#include "simulation.h"

#include <hpx/hpx.hpp>                  // hpx
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>

int hpx_main()
{
  hpx::cout << "Running VORTEX application - with HPX\n";

  double theta = 0.5;
  bool verify = true;
  double tol_ = 1e-8; // is supposed to be the same as tol in run_tests.cpp, but I was having linker errors

  double extT, mrtT, srtT, reoT, bldT, evaT, potT;
  int n, nnodes, NDST;
  int numtests = 3;
  bool printeach = true;

  run_test(extT, mrtT, srtT, reoT, bldT, evaT, potT,
           n, nnodes, NDST, theta, tol_, verify,
           numtests, printeach);
  return hpx::finalize();

}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);  // pass along command line arguments
}
