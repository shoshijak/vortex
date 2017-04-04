
#include <hpx/hpx.hpp>                  // hpx
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>

#include <unistd.h>                     // C
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>

#include "simulation.h"


int hpx_main()
{
    hpx::cout << "Starting HPX application - VORTEX with HPX\n";

    double theta = 0.5;
    bool verify = true;
    double tol_ = 1e-8; // is supposed to be the same as tol in run_tests.cpp, but I was having linker errors

    char filename[256];
  	strcpy(filename, "test-data/dN400");

  	if (access(filename, R_OK) == -1)
  	{
  		printf("WARNING: reference file <%s> not found.\n", filename);
  		return 1;
  	}
  	else
  		printf("reading from <%s> ...\n", filename);

  	FILE * fin = fopen(filename, "r");

  	assert(fin && sizeof(double) == sizeof(double));

  	test(theta, tol_, fin, verify);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);  // pass along command line arguments
}
