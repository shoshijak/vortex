#include "tree.h"

#define LMAX 15

void check(const double * ref, const double * res, const int N);

void evaluate(const Node* nodes, const double* expansions, const double *xdata, const double *ydata, const double *mdata,
		const double thetasquared, double * const result, const double xt, const double yt);

void potential(double theta, double *xsrc, double *ysrc, double *sources, int NSRC,
    double *xdst, double *ydst, int NDST, double *xtargets,
               double const& extT, double const& mrtT, double const& srtT,
               double const& reoT, double const& bldT, double const& evaT);

void test(double theta, double tol, FILE * f = NULL, bool verify = true, const double &extT, const double &mrtT, const double &srtT, const double &reoT, const double &bldT, const double &evaT, const double &potT);
