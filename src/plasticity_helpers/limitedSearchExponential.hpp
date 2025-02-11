#ifndef LIMITEDSEARCHEXPONENTIAL_HPP
#define LIMITEDSEARCHEXPONENTIAL_HPP

#include "../tools.hpp"


T limited_search_exponential(T epv, T epv_t, T epv_delta, T xi, T K, T p_t, T p00, T beta, int& right)
/**
 **** WRITTEN BY TOBIAS VERHEIJEN ****
 *
 * @param epv: Current plastic volumetric Hencky strain value
 * @param epv_t: Initial/trial plastic volumetric Hencky strain value
 * @param epv_delta: Full Newton-Raphson step direction
 * @param xi: xi value for simulation
 * @param K: Bulk modulus, (lambda + 2*mu/dim)
 * @param p_t: trial p: -K * Tr(epsilon)
 * @param p00: Initial consolidation pressure value, p0
 * @param beta: Cohesion parameter, indicates if the material can undergo tension without deformation (yield surface to the left of p=0)
 * @param right: Integer value indicating if the trial point started to the right or left of the apex of the yield surface
 * @return Updated epv value or infty if no valid step
 */
{
    T scale = 1.;
    int tries = 0;
    int max_tries = 32;
    int f = 0;
    T test_epv;

    // Ensures we stay within admissible region, this checks the next delta gamma value and ensures that it will be non-negative at the next step.
    while (tries < max_tries) {
        test_epv = epv - scale*epv_delta;  // Check epv value
        T pp0 = p00 * std::exp(-xi * test_epv);  // get next relative p0 value
        T pp = K*(test_epv - epv_t) + p_t;
        T val = 2. * pp + (beta - 1.) * pp0;  // compute denominator of delta gamma, only used for sign, removed K, Msq scaling terms to simplify computation
        T midpoint = (pp0 - beta * pp0) * 0.5;  // Find p-center of yield surface to avoid div-by-0 error

        // Check for compaction and dilation cases for correct direction of epv update.
        if (std::abs(pp - midpoint) <= T(1.e-5)) {
            break;
        }

        if (((right == 1 && (pp < midpoint || p_t < midpoint)) || (right == 0 && (pp > midpoint || p_t > midpoint))) && tries > -1) {
            test_epv = epv;
            right = -1;
            break;
        }

        if (sgn(epv_t - test_epv) == sgn(val) || sgn(epv_t - test_epv) == 0) {
            break;
        }
        scale *= 0.5;

        tries++;
        if (tries != max_tries) {
            continue;
        }

        // This indicates that no step in the original direction is admissible, check to see if step in other direction is possible
        tries = 0;
        scale = 1.;
        epv_delta *= -1.;  // invert direction
        f++;  // update flag

        if (f > 1){
            // indicates that no valid step in either direction can be made. Exit program
            return std::numeric_limits<double>::infinity();
        }
    }
    return test_epv;
}

#endif //LIMITEDSEARCHEXPONENTIAL_HPP
