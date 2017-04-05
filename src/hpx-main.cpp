
#include <hpx/hpx.hpp>                  // hpx
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>

#include "simulation.h"


int hpx_main()
{
#ifdef RUN_WITH_OMP
  hpx::cout << "Starting VORTEX application - with OpenMP\n";
#else
  hpx::cout << "Running VORTEX application - with HPX\n";
#endif

    double theta = 0.5;
    bool verify = true;
    double tol_ = 1e-8; // is supposed to be the same as tol in run_tests.cpp, but I was having linker errors

        double extT, mrtT, srtT, reoT, bldT, evaT, potT;
        int n, nnodes, NDST;
        int numtests = 50;
        bool printeach = false;

        run_test(extT, mrtT, srtT, reoT, bldT, evaT, potT,
                 n, nnodes, NDST, theta, tol_, verify,
                 numtests, printeach);


    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);  // pass along command line arguments
}
