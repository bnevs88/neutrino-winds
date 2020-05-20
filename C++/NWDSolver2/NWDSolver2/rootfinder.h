#ifndef FIND1DROOT_HPP_
#define FIND1DROOT_HPP_

#include <math.h>

// Searches for a root in the function func(x) between x = xlo and x = xup
// using a very simple bisection method 
template <class FUNCTION>
double Find1DRoot(double xlo, double xup, FUNCTION func, double tol = 1.e-8,
    bool verbose = false) {
    double fup = func(xup);
    double flo =0./*= func(xlo)*/;
    double xmid, fmid;
    if (verbose) std::cout << fup << " " << flo << std::endl;
    int iter = 0;
    do {
        ++iter;
        xmid = 0.5 * (xup + xlo);
        fmid = func(xmid);
        if (fmid * fup > 0.0) {
            xup = xmid;
            fup = fmid;
        }
        else {
            xlo = xmid;
            flo = fmid;
        }
        if (verbose) std::cout << xlo << " " << xmid << " " << fup << " " << flo << std::endl;
        if (iter > 100) {
            std::cout << "Exceeded max iteration count in rootfinder" << std::endl;
            return xmid;            
            //throw 0;
        }
    } while (fabs(fmid) > tol && fabs(xup / xlo - 1.0) > tol);
    if (verbose) std::cout << fabs(fmid) << ", " << fabs(xup / xlo - 1.0) << std::endl;
    return xmid;
}

#endif // FIND1DROOT_HPP_