#include "plastic_models.hpp"

bool ModifiedCamClayRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K)
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
            T r1 = pt - p - K    * delta_gamma * k;
            T r2 = qt - q - 3*mu * delta_gamma * l;

            if ( iter > 4 && std::abs(y) < 1e-3 && std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                // debug("ModifiedCamClayRMA: Breaking loop bc small rx at iter = ", iter);
                break;
            }
            if (iter == max_iter - 1){ // did not break loop
                debug("ModifiedCamClayRMA: FATAL did not exit loop at iter = ", iter);
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

            T J11 = -1 - K*delta_gamma*dkdp;
            T J13 = -K*k;

            T J22 = -1 - 3*mu*delta_gamma*dldq;
            T J23 = -3*mu*l;

            T J31 = k;
            T J32 = l;

            T det = J11*(-J32*J23) + J13*(-J31*J22);

            if (abs(det) < T(1e-6)){
                debug("ModifiedCamClayRMA: Determinant of Jacobian too small: det = ", det);
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

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}





bool ModifiedCamClayHardRMA(T& p, T& q, int& exit, T M, T epv, T beta, T mu, T K, T xi)
{
    T p0 = std::max(T(1e-3), K*std::sinh(-xi*epv));
    T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);

    if (y > 0) {

        T delta_gamma = 0;
        T pt = p;
        T qt = q;
        T ddydqq = 4*beta+2;

        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            T delta_epv = (p-pt) / K;
            p0          = std::max(T(1e-3),    -K * std::sinh(xi*(epv+delta_epv)));
            T ddp0dpp   = std::max(T(0), -xi*xi/K * std::sinh(xi*(epv+delta_epv)));
            T dp0dp     = 0;
            if ( (epv+delta_epv) < 0 )
                dp0dp   = -xi * std::cosh(xi*(epv+delta_epv));

            y        = M*M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            T dydp   = M*M * ( 2*p + (beta-1)*(dp0dp*p+p0) + 2*beta*p0*dp0dp );
            T ddydpp = M*M * ( 2 + (beta-1)*(p*ddp0dpp + 2*dp0dp) + 2*beta*(dp0dp*dp0dp + p0*ddp0dpp) );
            T dydq   = 2*q * (2*beta+1);

            T r1 = pt - p - K    * delta_gamma * dydp;
            T r2 = qt - q - 3*mu * delta_gamma * dydq;

            if ( iter > 4 && std::abs(y) < 1e-3 && std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                // debug("ModifiedCamClayRMA: Breaking loop bc small rx at iter = ", iter);
                break;
            }
            if (iter == max_iter - 1){ // did not break loop
                if (p0 > 1.01e-3){
                    debug("ModifiedCamClayHardRMA: FATAL did not exit loop at iter = ", iter);
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
                else{ // p0 too small
                    p = 1e-15;
                    q = 1e-15;
                    break;
                }
            }

            T J11 = -1 - K*delta_gamma*ddydpp;
            T J13 = -K*dydp;

            T J22 = -1 - 3*mu*delta_gamma*ddydqq;
            T J23 = -3*mu*dydq;

            T J31 = dydp;
            T J32 = dydq;

            T det = J11*(-J32*J23) + J13*(-J31*J22);

            if (abs(det) < T(1e-6)){
                debug("ModifiedCamClayHardRMA: Determinant of Jacobian too small: det = ", det);
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

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}




bool PerzynaMCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc)
{
    T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);

    if (y > 0) {

        T pt = p;
        T qt = q;
        T Cp =    K*dt / (p0*p0 * perzyna_visc);
        T Cq = 3*mu*dt / (p0*p0 * perzyna_visc);
        T dkdp = 2*M*M;
        T dldq = 4*beta+2;

        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            if (iter == max_iter - 1){ // did not break loop
                debug("PerzynaMCCRMA: FATAL did not exit loop at iter = ", iter);
                exit = 1;
            }

            y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            T k = M*M * (beta*p0 + 2*p - p0);
            T l = 2*q * (2*beta+1);
            T n = std::sqrt(k*k/d + 3.0/2.0*l*l);
            T r1 = pt - p - Cp * y*k/n;
            T r2 = qt - q - Cq * y*l/n;

            if ( iter > 4 && std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                // debug("PerzynaMCCRMA: Breaking the loop due to small residual");
                break;
            }

            T dndp =             k*dkdp / (n*d);
            T dndq = (3.0/2.0) * l*dldq / n;

            T tmp1 = (k*n - y*dndp) / (n*n);
            T tmp2 = (l*n - y*dndq) / (n*n);

            T Ja = -1 - Cp * (tmp1 * k + y/n * dkdp);
            T Jb =    - Cp * tmp2 * k;
            T Jc =    - Cq * tmp1 * l;
            T Jd = -1 - Cq * (tmp2 * l + y/n * dldq);

            T det = Ja*Jd - Jb*Jc;

            if (abs(det) <= T(1e-6)){
                debug("Determinant of Jacobian too small: det = ", det);
                p -= 0.001*r1;
                q -= 0.001*r2;
            } else{
                p -= ( Jd*r1 - Jb*r2) / det;
                q -= (-Jc*r1 + Ja*r2) / det;
            }

            if (q < 1e-15){
                q = 1e-15;
            }

            // debug(iter, ":  r1   = ", r1);
            // debug(iter, ":  r2   = ", r2);

        } // end for loop

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}

bool PerzynaSinterMCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T epv, T S, T visc, T Sinf, T tc, T ec, T xi)
{
    //debug("S = ", S);
    T mu_sqrt6 = mu*std::sqrt(6);
    T pc = std::max( T(0), p0 * (1 + std::sinh(-xi*epv)) * (1+S) );
    T y = (M * M * (p - pc) * (p + beta * pc) + (1 + 2 * beta) * (q * q)) / (p0*p0);

    if (y > 0) {

        T pt = p;
        T qt = q;

        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            T delta_epv = -(pt-p) / K;

            // T S     = std::max(T(0.0), (dt/tc*Sinf - (qt-q)/(mu_sqrt6*ec) ) / (1.0+dt/tc));
            T S     = std::max( T(0.0), Sinf - (qt-q)/(mu_sqrt6 * ec * (1.0+dt/tc)) );
            T dSdq  = 1.0 / ( (1+dt/tc)*ec*mu_sqrt6 );
            T dSdq2 = 0;

            T f     = std::max(T(1e-3),    p0 * (1 + std::sinh(-xi*(epv+delta_epv))));
            T dfdp  = 0;
            T dfdp2 = 0;
            if ( (epv+delta_epv) < 0 )
                dfdp  =       -(xi/K) * p0 * std::cosh(-xi*(epv+delta_epv));
                dfdp2 = (xi*xi/(K*K)) * p0 * std::sinh(-xi*(epv+delta_epv));

            T pc      = f     * (1.0 + S);
            T dpcdp   = dfdp  * (1.0 + S);
            T dpcdp2  = dfdp2 * (1.0 + S);
            T dpcdq   = f     * dSdq;
            T dpcdq2  = f     * dSdq2;
            T dpcdpdq = dfdp  * dSdq;

            T y      = (M * M * (p - pc) * (p + beta * pc) + (1 + 2 * beta) * (q * q)) / (p0*p0);
            T dydp   = - M*M/(p0*p0) * (2*beta*pc*dpcdp + (1-beta)*(dpcdp*p+pc)-2*p);
            T dydq   = 2*(1+2*beta)/(p0*p0) * q - M*M / (p0*p0) * (2*beta*pc*dpcdq + (1-beta)*p*dpcdq);
            T n      = std::sqrt(dydp*dydp/d + 3.0/2.0*dydq*dydq);
            T r1     = pt - p -    K*dt/visc * y * dydp/n;
            T r2     = qt - q - 3*mu*dt/visc * y * dydq/n;

            if (iter == max_iter - 1){ // did not break loop
                debug("PerzynaSinterMCCRMA: FATAL did not exit loop at iter = ", iter);
                debug("S  = ", S);
                debug("pc = ", pc);
                debug("y  = ", y);
                debug("r1 = ", r1);
                debug("r2 = ", r2);
                exit = 1;
            }

            // debug(iter, ":  r1   = ", r1);
            // debug(iter, ":  r2   = ", r2);

            if ( iter > 4 && std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                // debug("PerzynaMCCRMA: Breaking the loop due to small residual");
                break;
            }

            T dydp2  =                      - M*M/(p0*p0) * ( 2*beta*( dpcdp*dpcdp + pc*dpcdp2  ) + (1-beta)*( dpcdp2*p + dpcdp+dpcdp) - 2 );
            T dydq2  = 2*(1+2*beta)/(p0*p0) - M*M/(p0*p0) * ( 2*beta*( dpcdq*dpcdq + pc*dpcdq2  ) + (1-beta)*( p*dpcdq2 ) );
            T dydpdq =                      - M*M/(p0*p0) * ( 2*beta*( dpcdp*dpcdq + pc*dpcdpdq ) + (1-beta)*( dpcdq + p*dpcdpdq ) );

            T dndp = (2/d * dydp*dydp2  + 3 * dydq*dydpdq) / (2*n);
            T dndq = (2/d * dydp*dydpdq + 3 * dydq*dydq2) / (2*n);

            T dfracdp = (dydp * n - y * dndp) / (n*n);
            T dfracdq = (dydq * n - y * dndq) / (n*n);

            T Ja = -1 - K*dt/visc    * (dfracdp * dydp + y/n * dydp2);
            T Jb =    - K*dt/visc    * (dfracdq * dydp + y/n * dydpdq);
            T Jc =    - 3*mu*dt/visc * (dfracdp * dydq + y/n * dydpdq);
            T Jd = -1 - 3*mu*dt/visc * (dfracdq * dydq + y/n * dydq2);

            T det = Ja*Jd - Jb*Jc;

            if (abs(det) <= T(1e-6)){
                debug("Determinant of Jacobian too small: det = ", det);
                p -= 0.001*r1;
                q -= 0.001*r2;
            } else{
              // debug("Determinant acceptable: det = ", det);
                p -= ( Jd*r1 - Jb*r2) / det;
                q -= (-Jc*r1 + Ja*r2) / det;
            }

            if (q < 1e-15){
                q = 1e-15;
            }

        } // end for loop

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}
