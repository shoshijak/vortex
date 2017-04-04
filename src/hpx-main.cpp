
#include <hpx/hpx.hpp>                  // hpx
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>

#include "simulation.h"


int hpx_main()
{
    hpx::cout << "Starting HPX application - VORTEX with HPX\n";

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

//    size_t N(25);
//    for(size_t i(0);i<N;i++)
      test(theta, tol_, fin, verify);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);  // pass along command line arguments
}


/*******************************************************************************

                  RESULTS

----------------- Regular code from TA (branch master)



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









*******************************************************************************/
