#include "plastic_models.hpp"

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





bool MCCHardRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T rma_prefac, T epv)
{
    T p0;

    //// ALT 1
    if (epv < 0){
        p0 = p00 * (1.0 - std::sinh(xi*epv));
    } else{
        p0 = p00 * (1.0 - std::tanh(xi*epv));
    }

    //// ALT 2
    // T shift_epv = -std::asinh(p00/K) / xi;
    // if (epv+shift_epv < 0){
    //     p0 = K*std::sinh(-xi*(epv+shift_epv));
    // } else{
    //     p0 = K*std::tanh(-xi*(epv+shift_epv));
    // }

    // p0 = std::max(T(1e-2), p0);

    T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
    // T y = M * M * (p - p0) * (p + beta * p0) + (q * q);

    if (y > 0) {

        T delta_gamma = 0;
        T pt = p;
        T qt = q;
        T ddydqq = 4*beta+2;
        // T ddydqq = 2;

        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            T delta_epv = (p-pt) / K;
            T dp0dp;
            T ddp0dpp;

            //// ALT 1
            T s1 = std::sinh(xi*(epv+delta_epv));
            T c1 = std::cosh(xi*(epv+delta_epv));
            T c2 = c1 * c1;
            T c3 = c1 * c1 * c1;
            if ((epv+delta_epv) < 0){
                p0      =                  p00 * (1.0 - s1);
                dp0dp   = -(xi/K)        * p00 * c1;
                ddp0dpp = -(xi*xi/(K*K)) * p00 * s1;
            } else{
                p0      =                  p00 * (1.0 - s1/c1);
                dp0dp   = -(xi/K)        * p00 / c2;
                ddp0dpp =  (xi*xi/(K*K)) * p00 * 2.0 * s1 / c3;
            }

            //// ALT 2
            // T s1 = std::sinh(-xi*(epv+shift_epv+delta_epv));
            // T c1 = std::cosh(-xi*(epv+shift_epv+delta_epv));
            // T c2 = c1 * c1;
            // T c3 = c1 * c1 * c1;
            // if ((epv+shift_epv+delta_epv) < 0){
            //     p0      =  K       * s1;
            //     dp0dp   = -xi      * c1;
            //     ddp0dpp =  xi*xi/K * s1;
            // } else{
            //     p0      =  K  * s1/c1;
            //     dp0dp   = -xi / c2;
            //     ddp0dpp = -xi*xi/K * 2.0*s1/c3;
            // }

            y        = M*M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            // y        = M*M * (p - p0) * (p + beta * p0) + (q * q);
            T dydp   = M*M * ( 2*p + (beta-1)*(dp0dp*p+p0) + 2*beta*p0*dp0dp );
            T ddydpp = M*M * ( 2 + (beta-1)*(p*ddp0dpp + 2*dp0dp) + 2*beta*(dp0dp*dp0dp + p0*ddp0dpp) );
            T dydq   = 2*q * (2*beta+1);
            // T dydq   = 2*q;

            T r1 = pt - p - K             * delta_gamma * dydp;
            T r2 = qt - q - rma_prefac*mu * delta_gamma * dydq;

            if ( iter > 4 && std::abs(y) < 1e-3 && std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                break;
            }
            if (iter == max_iter - 1){ // did not break loop
                if (p0 > 1.01e-1){
                    debug("RMA: FATAL did not exit loop at iter = ", iter);
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

            T J11 = -1 - K            *delta_gamma*ddydpp;
            T J13 = -K*dydp;

            T J22 = -1 - rma_prefac*mu*delta_gamma*ddydqq;
            T J23 = -rma_prefac*mu*dydq;

            T J31 = dydp;
            T J32 = dydq;

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

        } // end for loop

        p = std::max(p, -beta * p0);
        p = std::min(p, p0);
        q = M * std::sqrt((p0 - p) * (beta * p0 + p) / (1 + 2 * beta));
        // q = M * std::sqrt((p0 - p) * (beta * p0 + p));

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}

bool MCCHardExpRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T rma_prefac, T epv)
{
    T p0 = std::max(T(1e-2), p00 * std::exp(-xi*epv));
    // T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
    T y = M * M * (p - p0) * (p + beta * p0) + (q * q);

    if (y > 0) {

        T delta_gamma = 0;
        T pt = p;
        T qt = q;
        // T ddydqq = 4*beta+2;
        T ddydqq = 2;

        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            T delta_epv = (p-pt) / K;
            p0          = p00 * std::exp(-xi*(epv+delta_epv));
            T dp0dp     = -xi/K * p0;
            T ddp0dpp   = xi*xi / (K*K) * p0;

            // y        = M*M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            y        = M*M * (p - p0) * (p + beta * p0) + (q * q);
            T dydp   = M*M * ( 2*p + (beta-1)*(dp0dp*p+p0) + 2*beta*p0*dp0dp );
            T ddydpp = M*M * ( 2 + (beta-1)*(p*ddp0dpp + 2*dp0dp) + 2*beta*(dp0dp*dp0dp + p0*ddp0dpp) );
            // T dydq   = 2*q * (2*beta+1);
            T dydq   = 2*q;

            T r1 = pt - p - K             * delta_gamma * dydp;
            T r2 = qt - q - rma_prefac*mu * delta_gamma * dydq;

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

            T J11 = -1 - K            *delta_gamma*ddydpp;
            T J13 = -K*dydp;

            T J22 = -1 - rma_prefac*mu*delta_gamma*ddydqq;
            T J23 = -rma_prefac*mu*dydq;

            T J31 = dydp;
            T J32 = dydq;

            T det = J11*(-J32*J23) + J13*(-J31*J22);

            if (abs(det) < T(1e-6)){
                debug("RMA: Determinant of Jacobian too small: det = ", det, ", iter = ", iter);
                p           -= 0.001*r1;
                q           -= 0.001*r2;
                delta_gamma -= 0.001*y / (p00*p00);
            } else{
                p           -= ( -J32*J23*r1 + J32*J13*r2 - J22*J13*y ) / det;
                q           -= (  J31*J23*r1 - J31*J13*r2 - J11*J23*y ) / det;
                delta_gamma -= ( -J31*J22*r1 - J11*J32*r2 + J11*J22*y ) / det;
            }

            if (q < 1e-15){
                q = 1e-15;
            }

        } // end for loop

        p = std::max(p, -beta * p0);
        p = std::min(p, p0);
        // q = M * std::sqrt((p0 - p) * (beta * p0 + p) / (1 + 2 * beta));
        q = M * std::sqrt((p0 - p) * (beta * p0 + p));

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}






bool SinterMCCRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T dt, T Sc, T tc, T ec, T epv, T S)
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
