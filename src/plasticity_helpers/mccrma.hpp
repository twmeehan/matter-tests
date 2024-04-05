#ifndef MCCRMA_HPP
#define MCCRMA_HPP

#include "../tools.hpp"

bool MCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T rma_prefac)
{
    T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);

    if (y > 0) {

        T delta_gamma = 0;
        T pt = p;
        T qt = q;
        T dkdp = 2*M*M;
        T dldq = 4*beta+2;

        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            T k = M*M * (beta*p0 + 2*p - p0);
            T l = 2*q * (2*beta+1);
            T r1 = pt - p - K             * delta_gamma * k;
            T r2 = qt - q - rma_prefac*mu * delta_gamma * l;

            if ( iter > 4 && std::abs(y) < 1e-3 && std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                // debug("RMA: Breaking loop bc small rx at iter = ", iter);
                break;
            }
            if (iter == max_iter - 1){ // did not break loop
                debug("RMA: FATAL did not exit loop at iter = ", iter);
                debug(iter, ":  r1   = ", r1);
                debug(iter, ":  r2   = ", r2);
                debug(iter, ":  y    = ", y);
                debug(iter, ":  p0   = ", p0);
                debug(iter, ":  pt   = ", pt);
                debug(iter, ":  qt   = ", qt);
                debug(iter, ":  p    = ", p);
                debug(iter, ":  q    = ", q);
                exit = 1;
            }

            T J11 = -1 - K            *delta_gamma*dkdp;
            T J13 = -K*k;

            T J22 = -1 - rma_prefac*mu*delta_gamma*dldq;
            T J23 = -rma_prefac*mu*l;

            T J31 = k;
            T J32 = l;

            T det = J11*(-J32*J23) + J13*(-J31*J22);

            if (abs(det) < T(1e-6)){
                debug("RMA: Determinant of Jacobian too small: det = ", det, ", iter = ", iter);
                p           -= 0.001*r1;
                q           -= 0.001*r2;
                delta_gamma -= 0.001*y / (p0*p0);
            } else{
                p           -= ( -J32*J23*r1 + J32*J13*r2 - J22*J13*y ) / det;
                q           -= (  J31*J23*r1 - J31*J13*r2 - J11*J23*y ) / det;
                delta_gamma -= ( -J31*J22*r1 - J11*J32*r2 + J11*J22*y ) / det;
            }

            if (q < 1e-15){
                q = 1e-15;
            }
            // debug(iter, ":  r1   = ", r1);
            // debug(iter, ":  r2   = ", r2);
            // debug(iter, ":  y    = ", y);

        } // end for loop

        p = std::max(p, -beta * p0);
        p = std::min(p, p0);
        q = M * std::sqrt((p0 - p) * (beta * p0 + p) / (1 + 2 * beta));

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}


#endif // MCCRMA_HPP
