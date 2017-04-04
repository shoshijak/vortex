
#include <hpx/hpx.hpp>                  // hpx
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>

#include <unistd.h>                     // C
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>

#include <algorithm>                    // STL
#include <limits>
#include <vector>

#include "timer.h"                      // vrtx
#include "tree.h"
#include "kernels.h"

int hpx_main()
{
    hpx::cout << "Starting HPX application - VORTEX with HPX\n";

    //!

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);  // pass along command line arguments
}
