#include <unistd.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>

#include <algorithm>
#include <limits>
#include <vector>

#include <hpx/include/parallel_for_loop.hpp>

#include "timer.h"
#include "tree.h"
#include "kernels.h"
#include "simulation.h"

double tol = 1e-8;


void check(const double * ref, const double * res, const int N)
{
	double linf = 0, l1 = 0, linf_rel = 0, l1_rel = 0;

	for(int i = 0; i < N; ++i)
	{
		assert(!std::isnan(ref[i]));
		assert(!std::isnan(res[i]));

		const double err = ref[i] - res[i];
		const double maxval = std::max(fabs(res[i]), fabs(ref[i]));
		const double relerr = err/std::max(1e-6, maxval);

		if (fabs(relerr) >= tol && fabs(err) >= tol)
			printf("%d: %e ref: %e -> %e %e\n", i, res[i], ref[i], err, relerr);

		assert(fabs(relerr) < tol || fabs(err) < tol);

		l1 += fabs(err);
		l1_rel += fabs(relerr);

		linf = std::max(linf, fabs(err));
		linf_rel = std::max(linf_rel, fabs(relerr));
	}

	printf("l-infinity errors: %.03e (absolute) %.03e (relative)\n", linf, linf_rel);
	printf("       l-1 errors: %.03e (absolute) %.03e (relative)\n", l1, l1_rel);
}


void evaluate(const Node* nodes, const double* expansions, const double *xdata, const double *ydata, const double *mdata,
		const double thetasquared, double * const result, const double xt, const double yt)
{
	enum { BUFSIZE = 16 };

	int stack[LMAX * 3];

	int bufcount = 0;
	double rzs[BUFSIZE], izs[BUFSIZE], masses[BUFSIZE];
	const double *rxps[BUFSIZE], *ixps[BUFSIZE];

	int stackentry = 0, maxentry = 0;

	stack[0] = 0;
	*result = 0;
	while(stackentry > -1)
	{
		const int nodeid = stack[stackentry--];
		const Node * const node = nodes + nodeid;

		assert(nodeid < node->child_id || node->child_id == 0);

		const double r2 = (xt - node->xcom)*(xt - node->xcom) + (yt - node->ycom)*(yt - node->ycom);

		if (node->r * node->r < thetasquared * r2)
		{
			rzs[bufcount] = xt - node->xcom;
			izs[bufcount] = yt - node->ycom;
			masses[bufcount] = node->mass;
			rxps[bufcount] = expansions + ORDER * (2 * nodeid + 0);
			ixps[bufcount] = expansions + ORDER * (2 * nodeid + 1);
			++bufcount;

			if (bufcount == BUFSIZE)
			{
				bufcount = 0;
				*result += e2p(rzs, izs, masses, rxps, ixps, BUFSIZE);
			}
		}
		else
		{
			if (node->child_id == 0)
			{
				const int s = node->part_start;
				const int e = node->part_end;

				*result += p2p(xdata + s, ydata + s, mdata + s, e - s, xt, yt);
			}
			else
			{
				for(int c = 0; c < 4; ++c)
					stack[++stackentry] = node->child_id + c;

				maxentry = std::max(maxentry, stackentry);
			}
		}
	}

	if (bufcount)
		*result += e2p(rzs, izs, masses, rxps, ixps, bufcount);

}


void potential(double theta, double *xsrc, double *ysrc, double *sources, int NSRC,
			  double *xdst, double *ydst, int NDST, double *xtargets)
{
	double *x, *y, *m, *xsorted, *ysorted, *msorted;
	double *expansions;

	x = xsrc;
	y = ysrc;
	m = sources;
	int n = NSRC;

	posix_memalign((void **)&xsorted, 32, sizeof(double) * n);
	posix_memalign((void **)&ysorted, 32, sizeof(double) * n);
	posix_memalign((void **)&msorted, 32, sizeof(double) * n);

	int k = 32*1;	// leaf capacity
	int maxnodes = (n + k - 1) / k * 60;
    Node* nodes;
	posix_memalign((void **)&nodes, 32, sizeof(Node) * maxnodes);
	posix_memalign((void **)&expansions, 32, sizeof(double) * 2 * ORDER * maxnodes);

	build(x, y, m, n, k, xsorted, ysorted, msorted, nodes, expansions);

	double thetasquared = theta*theta;

	Timer tm;
	tm.start();

#ifdef RUN_WITH_OMP
	#pragma omp parallel for schedule(static,1)
	for(int i = 0; i < NDST; ++i)
	  {
		  evaluate(nodes, expansions, xsorted, ysorted, msorted,
				  thetasquared, xtargets + i, xdst[i], ydst[i]);
	  }
#else
	hpx::parallel::for_loop(
		hpx::parallel::execution::par, 0, NDST,
		[&](int i)
#endif
	{
		evaluate(nodes, expansions, xsorted, ysorted, msorted,
				thetasquared, xtargets + i, xdst[i], ydst[i]);
	});

	double t = tm.elapsed();
	printf("Evaluation took %.3f ms (%.3f us per target)\n", t*1e-6, t*1e-3 / NDST);

	free(xsorted);
	free(ysorted);
	free(msorted);
	free(nodes);
	free(expansions);
}

void test(double theta, double tol, FILE * f, bool verify)
{
	int NSRC;
	fread(&NSRC, sizeof(int), 1, f);

	double *xsrc, *ysrc, *sources;
	posix_memalign((void **)&xsrc, 32, sizeof(double) * NSRC);
	posix_memalign((void **)&ysrc, 32, sizeof(double) * NSRC);
	posix_memalign((void **)&sources, 32, sizeof(double) * NSRC);

	fread(xsrc, sizeof(double), NSRC, f);
	fread(ysrc, sizeof(double), NSRC, f);
	fread(sources, sizeof(double), NSRC, f);

	int NDST;
	fread(&NDST, sizeof(int), 1, f);

	double *xdst, *ydst, *xref, *yref;
	posix_memalign((void **)&xdst, 32, sizeof(double) * NDST);
	posix_memalign((void **)&ydst, 32, sizeof(double) * NDST);
	posix_memalign((void **)&xref, 32, sizeof(double) * NDST);
	posix_memalign((void **)&yref, 32, sizeof(double) * NDST);

	fread(xdst, sizeof(double), NDST, f);
	fread(ydst, sizeof(double), NDST, f);
	fread(xref, sizeof(double), NDST, f);

	const double eps = std::numeric_limits<double>::epsilon() * 10;

	double *xtargets, *ytargets;
	posix_memalign((void **)&xtargets, 32, sizeof(double) * NDST);
	posix_memalign((void **)&ytargets, 32, sizeof(double) * NDST);

	printf("Testing %s with %d sources and %d targets (theta %.3e)...\n", "POTENTIAL", NSRC, NDST, theta);

	Timer tm;
	tm.start();
	potential(theta, xsrc, ysrc, sources, NSRC, xdst, ydst, NDST, xtargets);
	const double tpot = tm.elapsed();

	printf("\x1b[94msolved in %.2f ms\x1b[0m\n", tpot * 1e-6);

	if (verify)
	{
		const int OFFSET = 0;
		const int JUMP = 1;
		for(int i = OFFSET; i < NDST; i += JUMP)
		{
			const double xd = xdst[i];
			const double yd = ydst[i];

			double s = 0;

			for(int j = 0; j < NSRC; ++j)
			{
				const double xr = xd - xsrc[j];
				const double yr = yd - ysrc[j];
				const double r2 = xr * xr + yr * yr;
				const double f  = fabs(r2) > eps;
				s += 0.5 * f * log(r2 + eps) * sources[j];
			}
			xref[i] = s;
		}

		std::vector<double> a, b, c, d;

		for(int i = OFFSET; i < NDST; i += JUMP)
		{
			a.push_back(xref[i]);
			b.push_back(xtargets[i]);
			c.push_back(yref[i]);
			d.push_back(ytargets[i]);
		}

		check(&a[0], &b[0], a.size());
	}

	free(xdst);
	free(ydst);

	free(xtargets);
	free(ytargets);

	free(xref);
	free(yref);

	free(xsrc);
	free(ysrc);
	free(sources);

	printf("TEST PASSED.\n");
}
