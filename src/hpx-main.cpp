
#include "simulation.h"
#include <hpx/hpx.hpp>                  // hpx
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/parallel/executors.hpp>
#include <hpx/include/async.hpp>
#include <hpx/lcos/future_wait.hpp>
#include <boost/format.hpp>


int hpx_main(boost::program_options::variables_map& vm)
{
  // extract cmd line arg
  std::uint64_t hpx_task_threshold_ = vm["hpx_task_threshold"].as<std::uint64_t>();
  std::uint64_t hpx_task_threshold = hpx_task_threshold_;

  hpx::cout << "Running VORTEX application - with HPX\n"
            << "Task creation threshold (#particles) = " << hpx_task_threshold << "\n";

  double theta = 0.5;
  bool verify = true;
  double tol_ = 1e-8; // is supposed to be the same as tol in run_tests.cpp, but I was having linker errors

  double extT, mrtT, srtT, reoT, bldT, evaT, potT;
  int n, nnodes, NDST;
  int numtests = 1;
  bool printeach = true;

  run_test(extT, mrtT, srtT, reoT, bldT, evaT, potT,
           n, nnodes, NDST, theta, tol_, verify,
           numtests, printeach, hpx_task_threshold);

  return hpx::finalize();
}

int main(int argc, char* argv[])
{
    boost::program_options::options_description
        desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    desc_commandline.add_options() //! how to use: --hpx_task_threshold 30
            ("hpx_task_threshold",
             boost::program_options::value<std::uint64_t>()->default_value(0),
            "threshold for task creation in terms of # particles in the node");

    return hpx::init(desc_commandline, argc, argv);  // pass along command line arguments
}
