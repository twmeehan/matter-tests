#ifndef MCCSINTERRMA_HPP
#define MCCSINTERRMA_HPP

#include "../tools.hpp"

bool MCCSinterRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T dt, T Sc, T tc, T ec, T epv, T S)
{
    T p0 = p00 * std::exp(-xi*epv) * (1.0+S);
    T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);

    if (S < 0)
        debug("RMA: Negative S = ", S);

    if (y > 0) {

        T pt = p;
        T qt = q;

        T delta_gamma = 0;
        T S_orig = S;
        T mu_sqrt6 = mu * std::sqrt(6.0);
        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            T delta_epv = -(pt-p) / K;
            T delta_eps =  (qt-q) / mu_sqrt6;

            T S     = (S_orig + dt/tc*Sc - delta_eps/ec) / (1.0+dt/tc);
            T dSdq  = 1.0 / ( mu_sqrt6*ec*(1+dt/tc) );
            T ddSdqq = 0;

            T f     = p00 * std::exp(-xi*epv);
            T dfdp  = -xi/K * f;
            T ddfdpp = xi*xi / (K*K) * f;

            T p0      = f     * (1.0 + S);
            T dp0dp   = dfdp  * (1.0 + S);
            T ddp0dpp = ddfdpp * (1.0 + S);
            T dp0dq   = f     * dSdq;
            T ddp0dqq = f     * ddSdqq;
            T dp0dpdq = dfdp  * dSdq;

            y        = M*M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            T dydp   = M*M * ( 2*p + (beta-1)*(dp0dp*p+p0)          + 2*beta*p0*dp0dp );
            T ddydpp = M*M * ( 2   + (beta-1)*(p*ddp0dpp + 2*dp0dp) + 2*beta*(dp0dp*dp0dp + p0*ddp0dpp) );
            T dydq   = 2*(1+2*beta)*q + M*M * ( (beta-1)*p*dp0dq   - 2*beta*p0*dp0dq );
            T ddydqq = 2*(1+2*beta)   + M*M * ( (beta-1)*p*ddp0dqq - 2*beta*(dp0dq*dp0dq + p0*ddp0dqq) );
            T ddydpdq = 0; // to be changed!

            T r1 = pt - p - K    * delta_gamma * dydp;
            T r2 = qt - q - 3*mu * delta_gamma * dydq;

            if ( iter > 4 && std::abs(y) < 1e-3 && std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                break;
            }
            if (iter == max_iter - 1){ // did not break loop
                if (p0 > 1.01e-1){
                    debug("RMA: FATAL did not exit loop at iter = ", iter, ", iter = ", iter);
                    debug(iter, ":  r1   = ", r1);
                    debug(iter, ":  r2   = ", r2);
                    debug(iter, ":  y    = ", y);
                    debug(iter, ":  p0   = ", p0);
                    debug(iter, ":  pt   = ", pt);
                    debug(iter, ":  qt   = ", qt);
                    debug(iter, ":  p    = ", p);
                    debug(iter, ":  q    = ", q);
                    // exit = 1;
                }
                else{ // p0 too small
                    p = 1e-15;
                    q = 1e-15;
                    break;
                }
            }

            T J11 = -1 - K*delta_gamma*ddydpp;
            T J12 = -K*delta_gamma*ddydpdq;
            T J13 = -K*dydp;

            T J21 = -3*mu*delta_gamma*ddydpdq;
            T J22 = -1 - 3*mu*delta_gamma*ddydqq;
            T J23 = -3*mu*dydq;

            T J31 = dydp;
            T J32 = dydq;
            T J33 = 0;

            T det = J11*(-J32*J23) + J13*(-J31*J22); // to be changed

            if (abs(det) < T(1e-6)){
                debug("RMA: Determinant of Jacobian too small: det = ", det, ", iter = ", iter);
                p           -= 0.001*r1;
                q           -= 0.001*r2;
                delta_gamma -= 0.001*y / (p00*p00);
            } else{
                p           -= ( -J32*J23*r1 + J32*J13*r2 - J22*J13*y ) / det; // to be changed
                q           -= (  J31*J23*r1 - J31*J13*r2 - J11*J23*y ) / det; // to be changed
                delta_gamma -= ( -J31*J22*r1 - J11*J32*r2 + J11*J22*y ) / det; // to be changed
            }

            if (q < 1e-15){
                q = 1e-15;
            }

        } // end for loop

        p = std::max(p, -beta * p0);
        p = std::min(p, p0);
        q = M * std::sqrt((p0 - p) * (beta * p0 + p) / (1 + 2 * beta));

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}

#endif // MCCSINTERRMA_HPP
