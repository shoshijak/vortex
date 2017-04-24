#include <omp.h>

#include <cstdio>

#if  !defined(__INTEL_COMPILER) && !defined(__clang__)
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

#ifndef RUN_WITH_OMP
#include <hpx/parallel/algorithms/minmax.hpp>
#include <hpx/include/parallel_for_loop.hpp>
#include <hpx/parallel/algorithms/sort.hpp>
#include <hpx/parallel/executor_parameters.hpp>
#include <hpx/include/lcos.hpp>
#endif

#include <cmath>

#define LMAX 15 // probs could do much smaller - also defined in simulation.h

const double EPS = 10000 * std::numeric_limits<double>::epsilon() ;

void minmax_vec(const double xsrc[], const double ysrc[], const int nsources, double xmin_xmax_ymin_ymax[])
{
	double lxmi = 1e13, lymi = 1e13, lxma = -1e13, lyma = -1e13;

#ifdef RUN_WITH_OMP
	for(int i = 0; i < nsources; ++i)
	{
		const double xval = xsrc[i];
		const double yval = ysrc[i];

		lxmi = fmin(lxmi, xval);
		lxma = fmax(lxma, xval);

		lymi = fmin(lymi, yval);
		lyma = fmax(lyma, yval);
	}

        xmin_xmax_ymin_ymax[0] = lxmi;
        xmin_xmax_ymin_ymax[1] = lxma;
        xmin_xmax_ymin_ymax[2] = lymi;
        xmin_xmax_ymin_ymax[3] = lyma;
#else
        auto minmaxX = hpx::parallel::minmax_element(
                    hpx::parallel::execution::par,
                    xsrc,
                    xsrc+nsources);
        auto minmaxY = hpx::parallel::minmax_element(
                    hpx::parallel::execution::par,
                    ysrc,
                    ysrc+nsources);

        xmin_xmax_ymin_ymax[0] = *minmaxX.first;
        xmin_xmax_ymin_ymax[1] = *minmaxX.second;
        xmin_xmax_ymin_ymax[2] = *minmaxY.first;
        xmin_xmax_ymin_ymax[3] = *minmaxY.second;
#endif

}

void extent(const int N, const double* const x, const double* const y,
			double& xmin, double& ymin, double& ext)
{
#ifdef RUN_WITH_OMP

        static const int chunksize = 1024 * 4;
        double ext0, ext1;
	
	{
		const int nthreads = omp_get_max_threads();

		double xpartials[2][nthreads], ypartials[2][nthreads];

#pragma omp parallel
		{
			const int tid = omp_get_thread_num();

			double lxmi = 1e13, lymi = 1e13, lxma = -1e13, lyma = -1e13;

			const int start = tid * chunksize;
			const int step = nthreads * chunksize;

			for(int i = start; i < N; i += step)
			{
				double result[4];

				minmax_vec(x + i, y + i, std::min((int)chunksize, N - i), result);

				lxmi = std::min(lxmi, result[0]);
				lxma = std::max(lxma, result[1]);
				lymi = std::min(lymi, result[2]);
				lyma = std::max(lyma, result[3]);
			}

			xpartials[0][tid] = lxmi;
			xpartials[1][tid] = lxma;
			ypartials[0][tid] = lymi;
			ypartials[1][tid] = lyma;
		}

		xmin = *std::min_element(xpartials[0], xpartials[0] + nthreads);
		ymin = *std::min_element(ypartials[0], ypartials[0] + nthreads);

        ext0 = (*std::max_element(xpartials[1], xpartials[1] + nthreads) - xmin);
        ext1 = (*std::max_element(ypartials[1], ypartials[1] + nthreads) - ymin);
	}
#else
#ifdef PRINT
    std::cout << "[tree_prepare.h] in function extent" << std::endl;
#endif
         hpx::parallel::static_chunk_size param;
         hpx::parallel::execution::parallel_task_policy par_policy;
         auto policy = par_policy.with(param);
#ifdef PRINT
    std::cout << "[tree_prepare.h] in function extent: about to launch Minmax element search." << std::endl;
#endif

         auto minmaxX_ = hpx::parallel::minmax_element(policy, x, x+N);
#ifdef PRINT
    std::cout << "[tree_prepare.h] in function extent: Minmax element search launched on X." << std::endl;
#endif

         auto minmaxY_ = hpx::parallel::minmax_element(
                    policy, y, y+N);
#ifdef PRINT
    std::cout << "[tree_prepare.h] in function extent: Minmax element search launched on Y." << std::endl;
#endif

                auto minmaxX = minmaxX_.get();
                auto minmaxY = minmaxY_.get();
                xmin = *minmaxX.first;
                ymin = *minmaxY.first;
#ifdef PRINT
    std::cout << "[tree_prepare.h] in function extent: Minmax elements found." << std::endl;
#endif
                double ext0 = *minmaxX.second - xmin;
                double ext1 = *minmaxY.second - ymin;
#endif

	// For numerical reasons, shift the domain boundaries out by epsilon
        ext = fmax(ext0, ext1) * (1 + 2 * EPS);
        xmin -= EPS * ext;
        ymin -= EPS * ext;
#ifdef PRINT
    std::cout << "[tree_prepare.h] in function extent: leaving." << std::endl;
#endif
}

void morton(const int N, const double* const x, const double* const y,
			const double xmin, const double ymin, const double ext, int* index)
{

#ifdef RUN_WITH_OMP
	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		int xid = floor((x[i] - xmin) / ext * (1 << LMAX));
		int yid = floor((y[i] - ymin) / ext * (1 << LMAX));

		xid = (xid | (xid << 8)) & 0x00FF00FF;
		xid = (xid | (xid << 4)) & 0x0F0F0F0F;
		xid = (xid | (xid << 2)) & 0x33333333;
		xid = (xid | (xid << 1)) & 0x55555555;

		yid = (yid | (yid << 8)) & 0x00FF00FF;
		yid = (yid | (yid << 4)) & 0x0F0F0F0F;
		yid = (yid | (yid << 2)) & 0x33333333;
		yid = (yid | (yid << 1)) & 0x55555555;

		index[i] = xid | (yid << 1);
	}
#else
    hpx::parallel::for_loop(hpx::parallel::execution::par, 0, N,
                       [&](int i){
            int xid = floor((x[i] - xmin) / ext * (1 << LMAX));
            xid = (xid | (xid << 8)) & 0x00FF00FF;
            xid = (xid | (xid << 4)) & 0x0F0F0F0F;
            xid = (xid | (xid << 2)) & 0x33333333;
            xid = (xid | (xid << 1)) & 0x55555555;

            int yid = floor((y[i] - ymin) / ext * (1 << LMAX));
            yid = (yid | (yid << 8)) & 0x00FF00FF;
            yid = (yid | (yid << 4)) & 0x0F0F0F0F;
            yid = (yid | (yid << 2)) & 0x33333333;
            yid = (yid | (yid << 1)) & 0x55555555;

            index[i] = xid | (yid << 1);
    });
#endif
}

void sort(const int N, int* index, int* keys)
{
#ifdef PRINT
    std::cout << "sort: starting ..." << std::endl;
#endif
	std::pair<int, int> * kv = NULL;
	posix_memalign((void **)&kv, 32, sizeof(*kv) * N);
#ifdef PRINT
    std::cout << "sort: memory setup ..." << std::endl;
#endif
#ifdef RUN_WITH_OMP
	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		kv[i].first = index[i];
                kv[i].second = i;
	}
#else
        hpx::parallel::static_chunk_size param;
        auto policy = hpx::parallel::execution::par.with(param);
        hpx::parallel::for_loop(policy, 0, N,
                                [&](int i)
        {
                kv[i].first = index[i];
                kv[i].second = i;
        }
                                );
#ifdef PRINT
    std::cout << "sort: pairs generated" << std::endl;
#endif
#endif
#ifdef RUN_WITH_OMP
#if  !defined(__INTEL_COMPILER) && !defined(__clang__)
	__gnu_parallel::sort(kv, kv + N);
#else
	std::sort(kv, kv + N);
#endif
#else
        hpx::parallel::sort(
                    policy,
                    kv, kv+N);
#endif
#ifdef PRINT
    std::cout << "sort: parallel sort done" << std::endl;
#endif

#ifdef RUN_WITH_OMP
	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		index[i] = kv[i].first;
		keys[i] = kv[i].second;
	}
#else
        hpx::parallel::for_loop(policy,
                                0, N, [&](int i)
        {
                index[i] = kv[i].first;
                keys[i] = kv[i].second;
        }
                                );
#endif
#ifdef PRINT
    std::cout << "sort: indices and keys written" << std::endl;
#endif
	free(kv);
#ifdef PRINT
    std::cout << "sort: memory of pairs freed" << std::endl;
#endif
}

void reorder(const int N, const int* const keys, const double* const x, const double* const y, const double* const m, double* xsorted, double* ysorted, double *msorted)
{

#ifdef RUN_WITH_OMP
	#pragma omp parallel for
	for(int i = 0; i < N; ++i)
	{
		const int entry = keys[i];
		
		xsorted[i] = x[entry];
		ysorted[i] = y[entry];
		msorted[i] = m[entry];
	}
#else
    hpx::parallel::static_chunk_size param;
    auto policy = hpx::parallel::execution::par.with(param);
    hpx::parallel::for_loop(policy,
                            0, N, [&](int i)
    {
            const int entry = keys[i];
            xsorted[i] = x[entry];
            ysorted[i] = y[entry];
            msorted[i] = m[entry];
    });
#endif
}

/* Setup the node parameters in Treebuilder::build_tree
 * in       xsources[]
 *          ysources[]
 *          msources[]
 *          nsources
 * out      mass
 *          xcom
 *          ycom
 *          radius
 */
inline void node_setup(const double xsources[], const double ysources[], const double msources[], const int nsources,
                double& mass, double& xcom, double& ycom, double& radius)
{
#ifdef PRINT
    std::cout << "--node setup start--" << std::endl
              << "x " << xsources[0] << ", n " << nsources << std::endl
              << "m " << mass <<", xc " << xcom << ", r " << radius << std::endl;
#endif
        mass = 0;
        double weight = 0, xsum = 0, ysum = 0;

        for(int i = 0; i < nsources; ++i)
        {
                const double x = xsources[i];
                const double y = ysources[i];
                const double m = msources[i];
                const double w = fabs(m);

                mass += m;
                weight += w;
                xsum += x * w;
                ysum += y * w;
        }

        xcom = weight ? xsum / weight : 0;
        ycom = weight ? ysum / weight : 0;

        double r2 = 0;
        for(int i = 0; i < nsources; ++i)
        {
                const double xr = xsources[i] - xcom;
                const double yr = ysources[i] - ycom;

                r2 = fmax(r2, xr * xr + yr * yr);
        }

        radius = sqrt(r2);

#ifdef PRINT
        std::cout << "--node setup end--" << std::endl;
#endif
}


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int lower_bound_vec(int s, int e, const int val, const int keys[])
{
	int c = e - s;

	if (keys[s] >= val)
		return s;

	if (keys[e - 1] < val)
		return e;

	while (c)
	{
		const int s0 = s, e0 = e;

		const double h = (e - s) * 1.f / 8;

                for (int programIndex = 0; programIndex < 8; ++programIndex)
		{
			//int candidate_s = s0, candidate_e = e0;
			const int i = MIN(e0 - 1, (int)(s0 + programIndex * h + 0.499999f));

			const bool isless = keys[i] < val;
			const int candidate_s = isless ? i : s0;
			const int candidate_e = isless ? e0 : i;

			s = MAX(s, candidate_s);
			e = MIN(e, candidate_e);
		}
	
		c = MIN(c / 8, e - s);
	}

	return s + 1;
}

int upper_bound_vec(int s, int e, const int val, const int keys[])
{
	int c = e - s;

	if (keys[s] > val)
		return s;

	if (keys[e - 1] <= val)
		return e;

	while (c)
	{
		const int s0 = s, e0 = e;

		const double h = (e - s) * 1.f / 8;

		for(int programIndex = 0; programIndex < 8; ++programIndex)
		{
			//int candidate_s = s0, candidate_e = e0;
			const int i = MIN(e0 - 1, (int)(s0 + programIndex * h + 0.499999f));
	
			const bool isless = keys[i] <= val;
			const int candidate_s = isless ? i : s0;
			const int candidate_e = isless ? e0 : i;

			s = MAX(s, candidate_s);
			e = MIN(e, candidate_e);
		}
	
		c = MIN(c / 8, e - s);
	}

	return s + 1;
}

