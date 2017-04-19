#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include <atomic>
#include <utility>
#include <vector>
#include <functional>

#include "tree.h"
#include "kernels.h"
#include "tree_prepare.h"
#include "timer.h"

#include <omp.h>

#ifndef RUN_WITH_OMP
#include <hpx/include/lcos.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/actions.hpp>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <utility>
#include <string>

#include <boost/format.hpp>

#endif


class TreeBuilder
{
	Node   *nodes;
	double *expansions;

	int    maxnodes;
	double *xdata, *ydata, *mdata;
	int    *keydata;
	int K;
	
public:
#ifdef RUN_WITH_OMP
	int currnnodes;
#else
	std::atomic<int> currnnodes; //! if we're running with HPX, make this variable atomic
#endif

	TreeBuilder(Node* nodes, int maxnodes,
				double *xdata, double *ydata, double *mdata, int *keydata,
				double ext, double xmin, double ymin, int K, double *expansions) :
		nodes(nodes), maxnodes(maxnodes),
		xdata(xdata), ydata(ydata), mdata(mdata), keydata(keydata),
		K(K), currnnodes(1), expansions(expansions)
	{};

#ifdef RUN_WITH_OMP
	void build_tree_omp(const int nodeid); // using openMP tasking
#else
	hpx::future<void> build_tree_hpx(const int nodeid); // spawning HPX tasks
#endif
};

#ifdef RUN_WITH_OMP
void TreeBuilder::build_tree_omp(const int nodeid)
{
    //! TODO take this line out
    //! is this a safe access to currnnodes??
    std::cout << "Curr # nodes : " << currnnodes << ", nodeid : " << nodeid << std::endl;

	Node * const node = nodes + nodeid;

	const int s = node->part_start;
	const int e = node->part_end;
	const int l = node->level;
	const int mId = node->morton_id;
	const bool leaf = e - s <= K || l + 1 > LMAX;

    std::cout << " n : " << node << std::endl
              << " s : " << s << std::endl
              << " e : " << e << std::endl
              << " l : " << l << std::endl
              << " m : " << mId << std::endl
              << " ? : " << leaf << std::endl;
s
	node_setup(xdata + s, ydata + s, mdata + s, e - s,
				node->mass, node->xcom, node->ycom, node->r);
	p2e(xdata + s, ydata + s, mdata + s, e - s, node->xcom, node->ycom,
			expansions + 2*nodeid*ORDER, expansions + 2*nodeid*ORDER + ORDER);

    if(leaf == true)
    {
        std::cout << "--leaf--" << std::endl;
    } else if (leaf == false) {
        std::cout << "--split--" << std::endl;

	    int childbase;
	    #pragma omp atomic capture
	    {
		    childbase = currnnodes; currnnodes += 4;
	    }
	    assert(nodeid < childbase);
	    if (childbase + 4 >= maxnodes)
		    printf("node %d, chbase %d, maxnodes %d\n", nodeid,
			   childbase, maxnodes);
	    assert(childbase + 4 < maxnodes);

	    node->child_id = childbase;


		for(int c = 0; c < 4; ++c)
		{
			const int shift = 2 * (LMAX - l - 1);

			const int key1 = mId | (c << shift);
			const int key2 = key1 + (1 << shift) - 1;

			const size_t indexmin = s + c == 0 ? s : lower_bound_vec(s, e, key1, keydata);
			const size_t indexsup = s + c == 3 ? e : upper_bound_vec(s, e, key2, keydata);

			const int chId = childbase + c;

			nodes[chId].setup(indexmin, indexsup,
					  l + 1, key1, nodeid);

#pragma omp task firstprivate(chId) if (indexsup - indexmin > 5e3 && c < 3)
			{
				build_tree_omp(chId);
			}
		}
	}
}
#else
//! returns a future<void>
//! the build_tree_hpx(0) (ie the one called on the root node) is the "future that rules them all"
//! each execution of build_tree_hpx returns in one of these 2 ways:
//! 1. if leaf == true, it returns a future<void> which will be ready when p2e is ready
//! 2. if leaf == false, it returns a future<void> which will be ready when its own p2e, and the build_tree_hpx of each one of its children are ready
hpx::future<void> TreeBuilder::build_tree_hpx(const int nodeid)
{
    std::cout << "Curr # nodes : " << currnnodes.load() << ", curr nodeid : " << nodeid << std::endl;

    Node * const node = nodes + nodeid;

	const int  s    = node->part_start;
	const int  e    = node->part_end;
	const int  l    = node->level;
	const int  mId  = node->morton_id;
	const bool leaf = e - s <= K || l + 1 > LMAX; // K is 32

    std::cout << " n : " << node << std::endl << " s : " << s    << std::endl  << " e : " << e    << std::endl
              << " l : " << l    << std::endl << " m : " << mId  << std::endl  << " ? : " << leaf << std::endl;

	hpx::future<void> setup = hpx::async(
	      [&,node](){
	    node_setup(xdata + s, ydata + s, mdata + s, e - s,
		       node->mass, node->xcom, node->ycom, node->r);
	  });

	hpx::future<void> p2e_fut = setup.then(
            [&,node](hpx::future<void> fv)
            {
	            p2e(xdata + s, ydata + s, mdata + s, e - s,
                    node->xcom, node->ycom,
                    expansions + 2*nodeid*ORDER,
                    expansions + 2*nodeid*ORDER + ORDER);
                std::cout << "Here after p2e" << std::endl;
	  });

	if (leaf) {

        std::cout << "--leaf--" << std::endl;

//	    Stuff that's not quite what I want:
//	    p2e_fut.get(); // return hpx::future<void>;  //! John's suggestion... hmm... but it doesn't compile ... :/
//	    return hpx::async(p2e_fut); //! doesn't make sense
//	      What I want:
//	      1. return a future<void> that is ready when p2e_fut is ready
//	      2. non blocking
//            p2e_fut.get();
//	    return hpx::make_ready_future();
        return p2e_fut;
	}

	std::vector<hpx::future<void>> exp_and_child_fut;
	exp_and_child_fut.reserve(5);
	exp_and_child_fut.push_back(std::move(p2e_fut));

    //! atomically fetch the value of the last index in the array of nodes, and increment this index by 4 positions
	int childbase;
	childbase = currnnodes.fetch_add(4);

	//! for debugging purposes
	std::cout << "Childbase index is at : " << childbase << std::endl;

	assert(nodeid < childbase);
	if (childbase + 4 >= maxnodes)
		printf("node %d, chbase %d, maxnodes %d\n", nodeid, childbase, maxnodes);
	assert(childbase + 4 < maxnodes);

	node->child_id = childbase;

    for(int c = 0; c < 4; ++c)
    {
        const int shift = 2 * (LMAX - l - 1);

        const int key1 = mId | (c << shift);
        const int key2 = key1 + (1 << shift) - 1;

        const size_t indexmin = s + c == 0 ? s : lower_bound_vec(s, e, key1, keydata); //! put for_loop in lower_bound code
        const size_t indexsup = s + c == 3 ? e : upper_bound_vec(s, e, key2, keydata); //! put for_loop in upper_bound code

        const int chId = childbase + c;
        std::cout << "launching build for ChildId : " << chId << std::endl;
        nodes[chId].setup(indexmin, indexsup, l + 1, key1, nodeid);

        exp_and_child_fut.push_back(
                hpx::async(
                        build_tree_hpx,chId
          ));
    }

    std::cout << "--children all launched--" << exp_and_child_fut.size() << std::endl;

    //! (Asynchronously) wait for all the futures in the vector of futures
    //hpx::when_all(exp_and_child_fut).get();
    //return std::move(exp_and_child_fut[0]);
    return std::move(hpx::when_all(exp_and_child_fut)); //! this would be the "correct" thing to do. Just trying the other one for debugging
    //return hpx::make_ready_future();
}
#endif


#ifdef RUN_WITH_OMP
// n = nsrc, k = max leaf capacity
void build(const double* const x, const double*const y, const double* mass, const int n, const int k,
	double* xsorted, double* ysorted, double* mass_sorted,
	   Node *nodes, double *expansions,
	   double& exTm, double& morTm, double& sorTm,
	   double& reoTm, double& bldTm, int& nnodes)
{
	Timer tm;
	int *index, *keys;
	double ext, xmin, ymin;

	posix_memalign((void **)&index, 32, sizeof(int) * n);
	posix_memalign((void **)&keys,  32, sizeof(int) * n);

	tm.start();
	extent(n, x, y, xmin, ymin, ext);
	exTm = 1e-6 * tm.elapsedAndReset();
	
	morton(n, x, y, xmin, ymin, ext, index);
	morTm = 1e-6 * tm.elapsedAndReset();
	
	sort(n, index, keys);
	sorTm = 1e-6 * tm.elapsedAndReset();
	
	reorder(n, keys, x, y, mass, xsorted, ysorted, mass_sorted);
	reoTm = 1e-6 * tm.elapsedAndReset();

#ifdef SMALL_SIMULATION
	if (n <= 100)
		for (int i=0; i<n; i++)
			printf("%2d : coordinates = [%6f %6f], ind = %x\n", i, xsorted[i], ysorted[i], index[i]); 
#endif
	auto builder = new TreeBuilder(nodes, (n + k - 1) / k * 6, //! printf?
						xsorted, ysorted, mass_sorted, index, ext, xmin, ymin, k, expansions);

	tm.start();
#ifdef RUN_WITH_OMP
	#pragma omp parallel
	#pragma omp single nowait
	{
		nodes[0].rootsetup(0, n, 0, 0);
		builder->build_tree_omp(0);
	}
	bldTm = 1e-6 * tm.elapsedAndReset();
	nnodes = builder->currnnodes;

	free(index);
	free(keys);
#else
	nodes[0].rootsetup(0, n, 0, 0);
	double dummyTm = 1e-6*tm.elapsedAndReset(); //! TODO put this here or one line lower??
	hpx::future<void> tree_fut = builder->build_tree_hpx(0);
	return tree_fut.then(hpx::launch::sync,
			     //! will only be called once the function build_tree_hpx(0) (run on the root node), which is effectively the "future that rules them all" is ready.
			     //! sync launch policy s.t. the timing is accurate: we call the lambda synchronously immediately after the future tree_fut is ready.
			     [&,builder](hpx::future<void> fv){ //! TODO review capture
	    bldTm = 1e-6 * tm.elapsedAndReset();
	    int a = builder->currnnodes;
	    nnodes = a;
            free(index);
            free(keys);
          });
#endif

}
#else
// n = nsrc, k = max leaf capacity
// returns a hpx::future<void>
// so that we move on to solving the system only once this future is ready
hpx::future<void> build(const double* const x, const double*const y, const double* mass, const int n, const int k,
                        double* xsorted, double* ysorted, double* mass_sorted,
                        Node *nodes, double *expansions,
                        double& exTm, double& morTm, double& sorTm,
                        double& reoTm, double& bldTm, int& nnodes)
{
    Timer tm;
    int *index, *keys;
    double ext, xmin, ymin;

    posix_memalign((void **)&index, 32, sizeof(int) * n);
    posix_memalign((void **)&keys,  32, sizeof(int) * n);

    tm.start();
    extent(n, x, y, xmin, ymin, ext);
    exTm = 1e-6 * tm.elapsedAndReset();

    morton(n, x, y, xmin, ymin, ext, index);
    morTm = 1e-6 * tm.elapsedAndReset();

    sort(n, index, keys);
    sorTm = 1e-6 * tm.elapsedAndReset();

    reorder(n, keys, x, y, mass, xsorted, ysorted, mass_sorted);
    reoTm = 1e-6 * tm.elapsedAndReset();

#ifdef SMALL_SIMULATION
    if (n <= 100)
		for (int i=0; i<n; i++)
			printf("%2d : coordinates = [%6f %6f], ind = %x\n", i, xsorted[i], ysorted[i], index[i]);
#endif
    std::cout << "tree max # nodes : " << (n + k - 1) / k * 6 << std::endl;
    auto builder = new TreeBuilder(nodes, (n + k - 1) / k * 6,
                                   xsorted, ysorted, mass_sorted, index, ext, xmin, ymin, k, expansions);

    tm.start();
#ifdef RUN_WITH_OMP
    #pragma omp parallel
	#pragma omp single nowait
	{
		nodes[0].rootsetup(0, n, 0, 0);
		builder->build_tree_omp(0);
	}
	bldTm = 1e-6 * tm.elapsedAndReset();
	nnodes = builder->currnnodes;

	free(index);
	free(keys);
#else
    std::cout << "--setting up root node--" << std::endl;
    nodes[0].rootsetup(0, n, 0, 0);
    double dummyTm = 1e-6*tm.elapsedAndReset(); //! TODO put this here or one line lower??
    std::cout << "--launching build from root node--" << std::endl;
    hpx::future<void> tree_fut = builder->build_tree_hpx(0);
    std::cout << "--build launched from root node--" << std::endl;
    std::cout << "--waiting for future-that-rules-them-all to be ready--" << std::endl;
    return tree_fut.then(hpx::launch::sync,
            //! will only be called once the function build_tree_hpx(0) (run on the root node), which is effectively the "future that rules them all" is ready.
            //! sync launch policy s.t. the timing is accurate: we call the lambda synchronously immediately after the future tree_fut is ready.
                         [&,builder](hpx::future<void> fv){ //! TODO review capture
                             std::cout << "--future that rules them all ready--" << std::endl;
                             std::cout << "--tree building done, setting time--" << std::endl;
                             bldTm = 1e-6 * tm.elapsedAndReset();
                             int a = builder->currnnodes;
                             nnodes = a;
                             free(index);
                             free(keys);
                             std::cout << "--time set and memory freed!--" << std::endl;
                         });
#endif

}
#endif