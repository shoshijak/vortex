
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

    char filename[256];
        strcpy(filename, "../../test-data/dN400");

  	if (access(filename, R_OK) == -1)
  	{
  		printf("WARNING: reference file <%s> not found.\n", filename);
  		return 1;
  	}
  	else
  		printf("reading from <%s> ...\n", filename);

  	FILE * fin = fopen(filename, "r");

  	assert(fin && sizeof(double) == sizeof(double));

        double const& extT, mrtT, srtT, reoT, bldT, evaT, potT;

        test(theta, tol_, fin, verify,
             extT, mrtT, srtT, reoT, bldT, evaT, potT);

        printf("Evaluation took %.3f ms (%.3f us per target)\n", t*1e-6, t*1e-3 / NDST);
        printf("\x1b[94msolved in %.2f ms\x1b[0m\n", potT * 1e-6);


    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);  // pass along command line arguments
}


/*******************************************************************************

                  RESULTS

----------------  Regular code in an HPX main (branch vrtxhpx)

reading from <test-data/dN400> ...
Testing POTENTIAL with 216064 sources and 354 targets (theta 5.000e-01)...
TIME for N = 216064 (18085 nodes)  is  280.00 ms
	extent:  16.00 ms
	morton:  32.00 ms
	sorting: 124.00 ms
	reordering:  32.00 ms
	building:  76.00 ms
Evaluation took 17.337 ms (48.975 us per target)
solved in 345.34 ms
l-infinity errors: 5.504e-09 (absolute) 1.713e-06 (relative)
       l-1 errors: 2.597e-07 (absolute) 4.258e-05 (relative)
TEST PASSED.
Starting HPX application - VORTEX with HPX

OpenMP spawns extra threads: there is contention


----------------- Regular code from FS16 (branch master)

reading from <test-data/dN400> ...
Testing POTENTIAL with 216064 sources and 354 targets (theta 5.000e-01)...
TIME for N = 216064 (18085 nodes)  is   38.27 ms
	extent:   1.22 ms
	morton:   2.15 ms
	sorting:   5.50 ms
	reordering:   1.07 ms
	building:  28.33 ms
Evaluation took 3.640 ms (10.282 us per target)
solved in 46.23 ms
l-infinity errors: 5.504e-09 (absolute) 1.713e-06 (relative)
       l-1 errors: 2.597e-07 (absolute) 4.258e-05 (relative)
TEST PASSED.

----------------------------------- replaced (for, minmax, sort)

Testing POTENTIAL with 216064 sources and 354 targets (theta 5.000e-01)...
TIME for N = 216064 (18085 nodes)  is  158.40 ms
	extent:   1.20 ms
	morton:   1.25 ms
	sorting:  95.17 ms
	reordering:   3.54 ms
	building:  57.24 ms
Evaluation took 2.236 ms (6.316 us per target)
solved in 161.62 ms
l-infinity errors: 5.504e-09 (absolute) 1.713e-06 (relative)
       l-1 errors: 2.597e-07 (absolute) 4.258e-05 (relative)
TEST PASSED.
Running VORTEX application - with HPX


*******************************************************************************/
