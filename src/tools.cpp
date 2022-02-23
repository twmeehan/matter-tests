#include "tools.hpp"

unsigned int load_array(std::vector<TV>& array, std::string file_name)
{
    std::ifstream file(file_name);

    std::string line;
    T value;
    unsigned int p = 0; // particle
    unsigned int j;     // component (x, y or z)

    if ( file.is_open() ) {
        while ( std::getline(file, line) ) {
            j = 0;
            std::stringstream line_stream(line);
            while ( line_stream >> value ) {
                // debug("p = ", p, " j = ", j, " value = ", value);
                array[p](j) = value;
                j++;
            }
            p++;
        }
    }
    else {
        std::cout << "Unable to open '"<< file_name << "'" << std::endl;
    }

    return p;
}


// Taken from: https://gist.github.com/lorenzoriano/5414671
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}


bool copy_file(std::string source, std::string destination){
    std::ifstream in(source, std::ios::binary);
    std::ofstream out(destination, std::ios::binary);
    out << in.rdbuf();
    return in && out;
}




bool PerzynaNewCamClayRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc)
{
    typedef Eigen::Matrix<T, 2, 1> TV2; // 3 dim vector regardless of dim of problem
    using std::sqrt;
    using std::pow;

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
                debug("PerzynaNewCamClayRMA: FATAL did not exit loop at iter = ", iter);
                exit = 1;
            }

            y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            T k = M*M * (beta*p0 + 2*p - p0);
            T l = 2*q * (2*beta+1);
            T n = sqrt(k*k/d + 3.0/2.0*l*l);
            T r1 = pt - p - Cp * y*k/n;
            T r2 = qt - q - Cq * y*l/n;

            if ( std::abs(r1) < 1e-3 && std::abs(r2) < 1e-3 ){
                // debug("PerzynaCamClayRMA: Breaking the loop due to small residual");
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

            T p_prev = p;
            T q_prev = q;

            if (abs(det) <= T(1e-6)){
                p -= 0.001*r1;
                q -= 0.001*r2;
            } else{
                p -= ( Jd*r1 - Jb*r2) / det;
                q -= (-Jc*r1 + Ja*r2) / det;
            }

            if (q < 1e-10){
                q = 1e-10;
            }

            // debug(iter, ":  pt   = ", pt);
            // debug(iter, ":  qt   = ", qt);
            // debug(iter, ":  p    = ", p);
            // debug(iter, ":  q    = ", q);
            // debug(iter, ":  r1   = ", r1);
            // debug(iter, ":  r2   = ", r2);


            if ( iter > 4 && std::abs(p-p_prev) < 1e-3 && std::abs(q-q_prev) < 1e-3 ){
                // debug("PerzynaCamClayRMA: Breaking the loop due to no change in p and q");
                break;
            }
        } // end for loop

        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}



bool PerzynaCamClayRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc)
{
    typedef Eigen::Matrix<T, 2, 1> TV2; // 3 dim vector regardless of dim of problem
    using std::sqrt;
    using std::pow;

    T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);

    if (y > 0) {

        T pt = p;
        T qt = q;
        T v = perzyna_visc;

        T p_prev, q_prev;

        ///// PRECOMPS //////////
        T M_sq = M*M;
        T M_fo = M_sq * M_sq;
        T M_si = M_sq * M_fo;

        T p0_sq = p0*p0;
        T p0_fo = p0_sq*p0_sq;

        T sqrt_d = std::sqrt(d);
        T pow_d = d*sqrt_d;

        T tbpo = (2*beta + 1);
        T tbpo_sq = tbpo * tbpo;
        T tbpo_cu = tbpo * tbpo_sq;
        /////////////////////////

        TV2 r;
        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            if (iter == max_iter - 1){ // did not break loop
                debug("PerzynaCamClayRMA: FATAL did not exit loop at iter = ", iter);
                exit = 1;
            }

            y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);
            T dydp = M*M * (beta*p0 + 2*p - p0);
            T dydq = 2*q * (2*beta+1);
            T normalization = p0*p0;
            T prefactor = (dt/v) * (y/normalization) / sqrt(dydp*dydp/d + 3.0/2.0*dydq*dydq);
            T r0 = pt - p -    K*prefactor*dydp;
            T r1 = qt - q - 3*mu*prefactor*dydq;

            r = TV2( r0 , r1 );

            // ORIGINAL
            // T neg_det = (6*K*pow(M, 4)*d*pow(dt, 2)*mu*pow(q, 2)*pow(2*beta + 1, 2)*(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(-2*pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) - 12.0*d*pow(q, 2)*pow(2*beta + 1, 2) + 6.0*d*(2*beta + 1)*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1)))*(-pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + pow(M, 2)*(2*pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(4*beta + 2)) - 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*pow(beta*p0 + 2*p - p0, 2) - (-2*K*pow(M, 6)*sqrt(d)*dt*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1))*pow(beta*p0 + 2*p - p0, 2) + K*pow(M, 2)*sqrt(d)*dt*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 3.0L/2.0L)*(2*pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + pow(q, 2)*(4*beta + 2)) + pow(p0, 2)*v*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 2))*(-36.0*pow(d, 3.0L/2.0L)*dt*mu*pow(q, 2)*pow(2*beta + 1, 3)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1)) + 6*sqrt(d)*dt*mu*(2*beta + 1)*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 3.0L/2.0L)*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(6*beta + 3)) + pow(p0, 2)*v*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 2)))/(pow(p0, 4)*pow(v, 2)*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 4));
            //
            // T step0 = (K*pow(M, 2)*sqrt(d)*dt*q*(2*beta + 1)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(6*sqrt(d)*dt*mu*q*(2*beta + 1)*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1)) + pow(p0, 2)*v*(q - qt)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2)))*(-2*pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) - 12.0*d*pow(q, 2)*pow(2*beta + 1, 2) + 6.0*d*(2*beta + 1)*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1)))*(beta*p0 + 2*p - p0) + (K*pow(M, 2)*sqrt(d)*dt*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1))*(beta*p0 + 2*p - p0) + pow(p0, 2)*v*(p - pt)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2)))*(-36.0*pow(d, 3.0L/2.0L)*dt*mu*pow(q, 2)*pow(2*beta + 1, 3)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1)) + 6*sqrt(d)*dt*mu*(2*beta + 1)*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 3.0L/2.0L)*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(6*beta + 3)) + pow(p0, 2)*v*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 2)))/(pow(p0, 4)*pow(v, 2)*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 5.0L/2.0L));
            //
            // T step1 = (6*pow(M, 2)*sqrt(d)*dt*mu*q*(2*beta + 1)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(K*pow(M, 2)*sqrt(d)*dt*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1))*(beta*p0 + 2*p - p0) + pow(p0, 2)*v*(p - pt)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2)))*(-pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + pow(M, 2)*(2*pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(4*beta + 2)) - 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(beta*p0 + 2*p - p0) + (6*sqrt(d)*dt*mu*q*(2*beta + 1)*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1)) + pow(p0, 2)*v*(q - qt)*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2)))*(-2*K*pow(M, 6)*sqrt(d)*dt*sqrt(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2))*(pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(q, 2)*(2*beta + 1))*pow(beta*p0 + 2*p - p0, 2) + K*pow(M, 2)*sqrt(d)*dt*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 3.0L/2.0L)*(2*pow(M, 2)*(p - p0)*(beta*p0 + p) + pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + pow(q, 2)*(4*beta + 2)) + pow(p0, 2)*v*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 2)))/(pow(p0, 4)*pow(v, 2)*pow(pow(M, 4)*pow(beta*p0 + 2*p - p0, 2) + 6.0*d*pow(q, 2)*pow(2*beta + 1, 2), 5.0L/2.0L));

            // OPTIMIZED
            T q_sq = q*q;
            T Mpop = M_sq*(p - p0)*(beta*p0 + p);
            T tmp0 = (beta*p0 + 2*p - p0);
            T tmp1 = tmp0 * tmp0;
            T idiot = M_fo*tmp1 + 6.0*d*q_sq*tbpo_sq;
            T pre1 = sqrt(idiot);
            T pre2 = idiot * idiot;
            T pre3 = pre2 * pre1;// pow(M_fo*tmp1 + 6.0*d*q_sq*tbpo_sq, 5.0L/2.0L);

            T neg_det = (6*K*M_fo*d*dt*dt*mu*q_sq*tbpo_sq*idiot*(-2*M_fo*tmp1 - 12.0*d*q_sq*tbpo_sq + 6.0*d*tbpo*(Mpop + q_sq*tbpo))*(-M_fo*tmp1 + M_sq*(2*Mpop + q_sq*(4*beta + 2)) - 6.0*d*q_sq*tbpo_sq)*tmp1 - (-2*K*M_si*sqrt_d*dt*pre1*(Mpop + q_sq*tbpo)*tmp1 + K*M_sq*sqrt_d*dt*pow(M_fo*tmp1 + 6.0*d*q_sq*tbpo_sq, 1.5)*(2*Mpop + M_sq*tmp1 + q_sq*(4*beta + 2)) + p0_sq*v*pre2)*(-36.0*pow_d*dt*mu*q_sq*tbpo_cu*pre1*(Mpop + q_sq*tbpo) + 6*sqrt_d*dt*mu*tbpo*pow(M_fo*tmp1 + 6.0*d*q_sq*tbpo_sq, 1.5)*(Mpop + q_sq*(6*beta + 3)) + p0_sq*v*pre2))/(p0_fo*v*v*pow(M_fo*tmp1 + 6.0*d*q_sq*tbpo_sq, 4));

            T step0 = (K*M_sq*sqrt_d*dt*q*tbpo*pre1*(6*sqrt_d*dt*mu*q*tbpo*(Mpop + q_sq*tbpo) + p0_sq*v*(q - qt)*pre1)*(-2*M_fo*tmp1 - 12.0*d*q_sq*tbpo_sq + 6.0*d*tbpo*(Mpop + q_sq*tbpo))*tmp0 + (K*M_sq*sqrt_d*dt*(Mpop + q_sq*tbpo)*tmp0 + p0_sq*v*(p - pt)*pre1)*(-36.0*pow_d*dt*mu*q_sq*tbpo_cu*pre1*(Mpop + q_sq*tbpo) + 6*sqrt_d*dt*mu*tbpo*pow(M_fo*tmp1 + 6.0*d*q_sq*tbpo_sq, 1.5)*(Mpop + q_sq*(6*beta + 3)) + p0_sq*v*pre2))/(p0_fo*v*v*pre3);

            T step1 = (6*M_sq*sqrt_d*dt*mu*q*tbpo*pre1*(K*M_sq*sqrt_d*dt*(Mpop + q_sq*tbpo)*tmp0 + p0_sq*v*(p - pt)*pre1)*(-M_fo*tmp1 + M_sq*(2*Mpop + q_sq*(4*beta + 2)) - 6.0*d*q_sq*tbpo_sq)*tmp0 + (6*sqrt_d*dt*mu*q*tbpo*(Mpop + q_sq*tbpo) + p0_sq*v*(q - qt)*pre1)*(-2*K*M_si*sqrt_d*dt*pre1*(Mpop + q_sq*tbpo)*tmp1 + K*M_sq*sqrt_d*dt*pow(M_fo*tmp1 + 6.0*d*q_sq*tbpo_sq, 1.5)*(2*Mpop + M_sq*tmp1 + q_sq*(4*beta + 2)) + p0_sq*v*pre2))/(p0_fo*v*v*pre3);

            TV2 step( step0 , step1 );

            if (abs(neg_det) <= T(1e-6))
                step = T(-0.001) * r;
            else
                step = step / neg_det;

            p_prev = p;
            q_prev = q;

            // debug(iter, ": (p_prev, q_prev) = (", p_prev, ", ", q_prev,")");

            p += step(0);
            q += step(1);

            if ( iter > 4 && std::abs(p-p_prev) < 1e-3 && std::abs(q-q_prev) < 1e-3 ){
                // debug("PerzynaCamClayRMA: Breaking the loop");
                break;
            }
        } // end for loop
        if (q < 1e-10){
            q = 1e-10;
        }
        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}




bool PerzynaQuadRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc)
{
    typedef Eigen::Matrix<T, 2, 1> TV2; // 3 dim vector regardless of dim of problem
    using std::sqrt;
    using std::pow;

    T y = q * (1 + 2 * beta) + 2 * M * (p + beta * p0) * (p - p0) / p0;

    if (y > 0) {

        T pt = p;
        T qt = q;
        T v = perzyna_visc;

        T p_prev, q_prev;

        TV2 r;
        T max_iter = 40;
        for (int iter = 0; iter < max_iter; iter++) {

            // if (iter == max_iter - 1){ // did not break loop
            //     debug("PerzynaQuadRMA: FATAL did not exit loop at iter = ", iter);
            //     exit = 1;
            // }

            T dydp = 2*M/p0 * (beta*p0 + 2*p - p0);
            T dydq = 2*beta+1;
            T normalization = p0;
            T prefactor = (dt/v) * (y/normalization) / sqrt(dydp*dydp/d + 3/2*dydq*dydq);
            T r0 = pt - p -    K*prefactor*dydp;
            T r1 = qt - q - 3*mu*prefactor*dydq;

            r = TV2( r0 , r1 );

            // ORIGINAL
            // T neg_det = -(12*K*pow(M, 2)*d*pow(dt, 2)*mu*pow(2*beta + 1, 2)*(4*M*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + 6.0*pow(beta, 2)*d*pow(p0, 2) + 6.0*beta*d*pow(p0, 2) + 1.5*d*pow(p0, 2)) - pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 3.0L/2.0L))*pow(beta*p0 + 2*p - p0, 2) + (3*sqrt(d)*dt*mu*pow(2*beta + 1, 2) + v*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5)))*(-16*K*pow(M, 3)*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*pow(beta*p0 + 2*p - p0, 2)*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + 6.0*pow(beta, 2)*d*pow(p0, 2) + 6.0*beta*d*pow(p0, 2) + 1.5*d*pow(p0, 2)) + 4*K*M*sqrt(d)*dt*pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 3.0L/2.0L)*(2*M*(p - p0)*(beta*p0 + p) + M*pow(beta*p0 + 2*p - p0, 2) + p0*q*(2*beta + 1)) + pow(p0, 2)*v*pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 2)))/(pow(p0, 2)*pow(v, 2)*pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 5.0L/2.0L));
            //
            // T step0 = (-2*K*M*sqrt(d)*dt*(2*beta + 1)*(3*sqrt(d)*dt*mu*(2*beta + 1)*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1)) + p0*v*(q - qt)*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5)))*(beta*p0 + 2*p - p0) + (3*sqrt(d)*dt*mu*pow(2*beta + 1, 2) + v*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5)))*(2*K*M*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*(beta*p0 + 2*p - p0) + pow(p0, 2)*v*(p - pt)*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5))))/(pow(p0, 2)*pow(v, 2)*(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5)));
            //
            // T step1 = (6*M*sqrt(d)*dt*mu*(2*beta + 1)*(4*M*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + 6.0*pow(beta, 2)*d*pow(p0, 2) + 6.0*beta*d*pow(p0, 2) + 1.5*d*pow(p0, 2)) - pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 3.0L/2.0L))*(2*K*M*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*(beta*p0 + 2*p - p0) + pow(p0, 2)*v*(p - pt)*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5)))*(beta*p0 + 2*p - p0) + (3*sqrt(d)*dt*mu*(2*beta + 1)*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1)) + p0*v*(q - qt)*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5)))*(-16*K*pow(M, 3)*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*pow(beta*p0 + 2*p - p0, 2)*sqrt(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + 6.0*pow(beta, 2)*d*pow(p0, 2) + 6.0*beta*d*pow(p0, 2) + 1.5*d*pow(p0, 2)) + 4*K*M*sqrt(d)*dt*pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 3.0L/2.0L)*(2*M*(p - p0)*(beta*p0 + p) + M*pow(beta*p0 + 2*p - p0, 2) + p0*q*(2*beta + 1)) + pow(p0, 2)*v*pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 2)))/(pow(p0, 3)*pow(v, 2)*pow(4*pow(M, 2)*pow(beta*p0 + 2*p - p0, 2) + d*pow(p0, 2)*(2*beta + 1)*(3.0*beta + 1.5), 5.0L/2.0L));

            // OPTIMIZED
            T tmp1 = pow(beta*p0 + 2*p - p0, 2);
            T tmp = 4*M*M*tmp1 + d*p0*p0*(2*beta + 1)*(3.0*beta + 1.5);
            T tmp2 = pow(tmp, 3.0/2.0);
            T tmp3 = pow(tmp, 5.0/2.0);
            T tmp4 = tmp*tmp;
            T tmp5 = pow(2*beta + 1, 2);

            T neg_det = -(12*K*M*M*d*dt*dt*mu*tmp5*(4*M*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*sqrt(4*M*M*tmp1 + 6.0*beta*beta*d*p0*p0 + 6.0*beta*d*p0*p0 + 1.5*d*p0*p0) - tmp2)*tmp1 + (3*sqrt(d)*dt*mu*tmp5 + v*sqrt(tmp))*(-16*K*M*M*M*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*tmp1*sqrt(4*M*M*tmp1 + 6.0*beta*beta*d*p0*p0 + 6.0*beta*d*p0*p0 + 1.5*d*p0*p0) + 4*K*M*sqrt(d)*dt*tmp2*(2*M*(p - p0)*(beta*p0 + p) + M*tmp1 + p0*q*(2*beta + 1)) + p0*p0*v*tmp4))/(p0*p0*v*v*tmp3);

            T step0 = (-2*K*M*sqrt(d)*dt*(2*beta + 1)*(3*sqrt(d)*dt*mu*(2*beta + 1)*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1)) + p0*v*(q - qt)*sqrt(tmp))*(beta*p0 + 2*p - p0) + (3*sqrt(d)*dt*mu*tmp5 + v*sqrt(tmp))*(2*K*M*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*(beta*p0 + 2*p - p0) + p0*p0*v*(p - pt)*sqrt(tmp)))/(p0*p0*v*v*(tmp));

            T step1 = (6*M*sqrt(d)*dt*mu*(2*beta + 1)*(4*M*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*sqrt(4*M*M*tmp1 + 6.0*beta*beta*d*p0*p0 + 6.0*beta*d*p0*p0 + 1.5*d*p0*p0) - tmp2)*(2*K*M*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*(beta*p0 + 2*p - p0) + p0*p0*v*(p - pt)*sqrt(tmp))*(beta*p0 + 2*p - p0) + (3*sqrt(d)*dt*mu*(2*beta + 1)*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1)) + p0*v*(q - qt)*sqrt(tmp))*(-16*K*M*M*M*sqrt(d)*dt*(2*M*(p - p0)*(beta*p0 + p) + p0*q*(2*beta + 1))*tmp1*sqrt(4*M*M*tmp1 + 6.0*beta*beta*d*p0*p0 + 6.0*beta*d*p0*p0 + 1.5*d*p0*p0) + 4*K*M*sqrt(d)*dt*tmp2*(2*M*(p - p0)*(beta*p0 + p) + M*tmp1 + p0*q*(2*beta + 1)) + p0*p0*v*tmp4))/(p0*p0*p0*v*v*tmp3);

            TV2 step( step0 , step1 );

            if (abs(neg_det) <= T(1e-6))
                step = T(-0.001) * r;
            else
                step = step / neg_det;

            p_prev = p;
            q_prev = q;

            p += step(0);
            q += step(1);

            if ( iter > 4 && std::abs(p-p_prev) < 1e-3 && std::abs(q-q_prev) < 1e-3 ){
                // debug("PerzynaQuadRMA: Breaking the loop");
                break;
            }
        } // end for loop
        if (q < 1e-10){
            q = 1e-10;
        }
        return true; // if plastic, i.e., y > 0
    } // end if outside
    return false;
}





bool CamClayRMA(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus)
{
    typedef Eigen::Matrix<T, 3, 1> TV3; // 3 dim vector regardless of dim of problem
    using std::abs;
    using std::max;
    using std::min;
    using std::sqrt;

    T y = M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q);

    if (y > 0) {
        T max_p = p0;
        T min_p = -beta * max_p;
        T e = 1 + 2 * beta;
        T max_q = T(0.5) * M * max_p * (1 + beta) * (1.0/sqrt(e));

        p = max(min(p, max_p), min_p);
        q = min(q, max_q);
        if (max_q < T(1e-10)) {
            // Too small to project properly
            return true;
        }
        T scale = max(-min_p, max(max_p, max_q));
        T scale_inverse = 1 / scale;

        p0 *= scale_inverse;
        trace_epsilon *= scale_inverse;
        norm_eps_hat *= scale_inverse;
        p *= scale_inverse;
        q *= scale_inverse;

        // Project
        T a = 1 / (3 * mu);
        T b = 1 / bulk_modulus;
        T c = trace_epsilon;
        T d = -sqrt(T(2) / 3) * norm_eps_hat;
        T f = M * M;

        T gamma = 0;

        TV3 r;
        for (int iter = 0; iter < 40; iter++) {
            T d1 = (p - p0);
            T d2 = (p + beta * p0);
            T A13 = f * (d1 + d2);
            T A22 = a + 2 * e * gamma;
            r = TV3(c + b * p + A13 * gamma,
                d + q * A22,
                f * d1 * d2 + e * (q * q));

            T A11 = b + 2 * f * gamma;
            T A23 = 2 * e * q;
            T neg_det = (A13 * A13 * A22 + A11 * A23 * A23);
            TV3 step(-(A23 * A23 * r(0)) + A13 * A23 * r(1) - A13 * A22 * r(2),
                A13 * A23 * r(0) - A13 * A13 * r(1) - A11 * A23 * r(2),
                -(A13 * A22 * r(0)) - A11 * A23 * r(1) + A11 * A22 * r(2));

            if (abs(neg_det) <= T(1e-6))
                step = T(-0.001) * r;
            else
                step = step / neg_det;
            p += step(0);
            q += step(1);
            gamma += step(2);
        }
        p0 *= scale;
        p *= scale;
        q *= scale;

        p = max(min(p, max_p), min_p);
        q = min(abs(q), max_q);
        assert((M * M * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * (q * q)) <= T(1e-3));
        assert(std::isfinite(p));
        assert(std::isfinite(q));
        return true;
    }
    return false;
}

bool QuadRMA(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus)
{
    typedef Eigen::Matrix<T, 3, 1> TV3; // 3 dim vector regardless of dim of problem
    using std::abs;
    using std::max;
    using std::min;
    using std::sqrt;

    T y = q * (1 + 2 * beta) + 2 * M * (p + beta * p0) * (p - p0) / p0;

    if (y > 0) {

        T max_p = p0;
        T min_p = -beta * max_p;
        T max_q = M / (2 * beta + 1) * T(0.5) * p0 * (1 + beta) * (1 + beta); // do not use precomps here!

        p = max(min(p, max_p), min_p);
        q = min(q, max_q);
        if (max_q < T(1e-10)) {
            // Too small to project properly
            return true;
        }
        T scale = max(-min_p, max(max_p, max_q));
        T scale_inverse = 1 / scale;

        p0 *= scale_inverse;
        trace_epsilon *= scale_inverse;
        norm_eps_hat *= scale_inverse;
        p *= scale_inverse;
        q *= scale_inverse;

        // Precomputations
        T sqrt_6_normesphat = std::sqrt(6) * norm_eps_hat;
        T one_third = 1.0 / 3.0;
        T twobetaplusone = 2 * beta + 1;
        T twobetaplusone_sq = twobetaplusone * twobetaplusone;
        T p0_sq = p0 * p0;
        T M_sq = M * M;
        T twoM = 2 * M;
        T twoMK = twoM * bulk_modulus;
        T fourMK = 4 * M * bulk_modulus;
        T beta_p0 = beta * p0;
        T inv_p0 = 1.0 / p0;
        T inv_mu = 1.0 / mu;
        T inv_K = 1.0 / bulk_modulus;
        T inv_Kmu_p0_sq = 1.0 / (bulk_modulus * mu * p0_sq);
        T Kp0treps = bulk_modulus * p0 * trace_epsilon;
        T four_one_third_Msq_K = 4.0 * one_third * M_sq * bulk_modulus;

        T gamma = 0;
        TV3 r;
        for (int iter = 0; iter < 40; iter++) {

            T tmp = beta_p0 + 2 * p - p0; // NB contains p
            T tmp_sq = tmp * tmp;

            r = TV3(
                twoM * gamma * tmp * inv_p0 + trace_epsilon + p * inv_K,
                gamma * twobetaplusone - one_third * sqrt_6_normesphat + one_third * q * inv_mu,
                twoM * (p - p0) * (beta_p0 + p) * inv_p0 + q * twobetaplusone);

            T neg_det = (four_one_third_Msq_K * tmp_sq + mu * p0 * twobetaplusone_sq * (fourMK * gamma + p0)) * inv_Kmu_p0_sq;

            TV3 step(
                one_third * (twoMK * p0 * twobetaplusone * (mu * (3 * gamma * twobetaplusone - sqrt_6_normesphat) + q) * tmp - twoMK * (twoM * (p - p0) * (beta_p0 + p) + p0 * q * twobetaplusone) * tmp - 3 * mu * p0 * twobetaplusone_sq * (twoMK * gamma * tmp + Kp0treps + p * p0)) * inv_Kmu_p0_sq,
                (-four_one_third_Msq_K * (mu * (3 * gamma * twobetaplusone - sqrt_6_normesphat) + q) * tmp_sq + mu * twobetaplusone * (twoM * tmp * (twoMK * gamma * tmp + Kp0treps + p * p0) - (fourMK * gamma + p0) * (twoM * (p - p0) * (beta_p0 + p) + p0 * q * twobetaplusone))) * inv_Kmu_p0_sq,
                one_third * (-twoM * tmp * (twoMK * gamma * tmp + Kp0treps + p * p0) - p0 * twobetaplusone * (mu * (3 * gamma * twobetaplusone - sqrt_6_normesphat) + q) * (fourMK * gamma + p0) + (fourMK * gamma + p0) * (twoM * (p - p0) * (beta_p0 + p) + p0 * q * twobetaplusone)) * inv_Kmu_p0_sq);

            if (abs(neg_det) <= T(1e-6))
                step = T(-0.001) * r;
            else
                step = step / neg_det;
            p += step(0);
            q += step(1);
            gamma += step(2);

            if (  (q * (1 + 2 * beta) + 2 * M * (p + beta * p0) * (p - p0) / p0) <= T(1e-3)  ){ // NB: Unit scaled yield surface, this is thus a relative error
                break; // exit for loop
            }
            if (iter == 39){
                debug("QuadraticReturnMapping: FATAL did not exit loop at iter = ", iter);
            }
        }
        p0 *= scale;
        p *= scale;
        q *= scale;

        p = max(min(p, max_p), min_p);
        q = min(abs(q), max_q);
        assert((q * (1 + 2 * beta) + 2 * M * (p + beta * p0) * (p - p0) / p0) <= T(1e-3)); // yield surface, do not use precomps here! NB this is the unscaled yield surface!
        assert(std::isfinite(p));
        assert(std::isfinite(q));
        return true;
    }
    return false;
}


 // DO NOT USE - WRONG PROJECTION
bool QuadAnalyticRMA(T& p, T& q, int& exit, T M, T p0, T beta)
{
    using std::abs;
    using std::cbrt;
    using std::max;
    using std::min;
    using std::pow;
    using std::sqrt;

    T y = q * (1 + 2 * beta) + 2 * M * (p + beta * p0) * (p - p0) / p0;

    if (y > 0) {
        T max_p = p0;
        T min_p = -beta * max_p;
        T max_q = M / (2 * beta + 1) * T(0.5) * p0 * (1 + beta) * (1 + beta);

        T scale = max(-min_p, max(max_p, max_q));
        //T scale = 1.0;
        T scale_inverse = 1 / scale;
        p0 *= scale_inverse;
        p *= scale_inverse;
        q *= scale_inverse;

        // Lagrangian multiplier solution (The real solution)
        // T lamb = -1.0L/3.0L*(-3*(4*pow(M, 2)*pow(beta, 2)*pow(p0, 2) + 8*pow(M, 2)*beta*pow(p0, 2) + 4*pow(M, 2)*pow(p0, 2) - 16*M*beta*p0*q - 8*M*p0*q + 4*pow(beta, 2)*pow(p0, 2) + 4*beta*pow(p0, 2) + pow(p0, 2))/(16*pow(M, 2)*pow(beta, 2) + 16*pow(M, 2)*beta + 4*pow(M, 2)) + pow(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0, 2)/pow(4*M*pow(beta, 2) + 4*M*beta + M, 2))/cbrt((1.0L/2.0L)*sqrt(-4*pow(-3*(4*pow(M, 2)*pow(beta, 2)*pow(p0, 2) + 8*pow(M, 2)*beta*pow(p0, 2) + 4*pow(M, 2)*pow(p0, 2) - 16*M*beta*p0*q - 8*M*p0*q + 4*pow(beta, 2)*pow(p0, 2) + 4*beta*pow(p0, 2) + pow(p0, 2))/(16*pow(M, 2)*pow(beta, 2) + 16*pow(M, 2)*beta + 4*pow(M, 2)) + pow(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0, 2)/pow(4*M*pow(beta, 2) + 4*M*beta + M, 2), 3) + pow(27*(-2*M*beta*pow(p0, 3) + 2*M*beta*pow(p0, 2)*p - 2*M*pow(p0, 2)*p + 2*M*p0*pow(p, 2) + 2*beta*pow(p0, 2)*q + pow(p0, 2)*q)/(8*pow(M, 2)*pow(beta, 2) + 8*pow(M, 2)*beta + 2*pow(M, 2)) - 9*(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0)*(4*pow(M, 2)*pow(beta, 2)*pow(p0, 2) + 8*pow(M, 2)*beta*pow(p0, 2) + 4*pow(M, 2)*pow(p0, 2) - 16*M*beta*p0*q - 8*M*p0*q + 4*pow(beta, 2)*pow(p0, 2) + 4*beta*pow(p0, 2) + pow(p0, 2))/((4*M*pow(beta, 2) + 4*M*beta + M)*(16*pow(M, 2)*pow(beta, 2) + 16*pow(M, 2)*beta + 4*pow(M, 2))) + 2*pow(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0, 3)/pow(4*M*pow(beta, 2) + 4*M*beta + M, 3), 2)) + (27.0L/2.0L)*(-2*M*beta*pow(p0, 3) + 2*M*beta*pow(p0, 2)*p - 2*M*pow(p0, 2)*p + 2*M*p0*pow(p, 2) + 2*beta*pow(p0, 2)*q + pow(p0, 2)*q)/(8*pow(M, 2)*pow(beta, 2) + 8*pow(M, 2)*beta + 2*pow(M, 2)) - 9.0L/2.0L*(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0)*(4*pow(M, 2)*pow(beta, 2)*pow(p0, 2) + 8*pow(M, 2)*beta*pow(p0, 2) + 4*pow(M, 2)*pow(p0, 2) - 16*M*beta*p0*q - 8*M*p0*q + 4*pow(beta, 2)*pow(p0, 2) + 4*beta*pow(p0, 2) + pow(p0, 2))/((4*M*pow(beta, 2) + 4*M*beta + M)*(16*pow(M, 2)*pow(beta, 2) + 16*pow(M, 2)*beta + 4*pow(M, 2))) + pow(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0, 3)/pow(4*M*pow(beta, 2) + 4*M*beta + M, 3)) - 1.0L/3.0L*cbrt((1.0L/2.0L)*sqrt(-4*pow(-3*(4*pow(M, 2)*pow(beta, 2)*pow(p0, 2) + 8*pow(M, 2)*beta*pow(p0, 2) + 4*pow(M, 2)*pow(p0, 2) - 16*M*beta*p0*q - 8*M*p0*q + 4*pow(beta, 2)*pow(p0, 2) + 4*beta*pow(p0, 2) + pow(p0, 2))/(16*pow(M, 2)*pow(beta, 2) + 16*pow(M, 2)*beta + 4*pow(M, 2)) + pow(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0, 2)/pow(4*M*pow(beta, 2) + 4*M*beta + M, 2), 3) + pow(27*(-2*M*beta*pow(p0, 3) + 2*M*beta*pow(p0, 2)*p - 2*M*pow(p0, 2)*p + 2*M*p0*pow(p, 2) + 2*beta*pow(p0, 2)*q + pow(p0, 2)*q)/(8*pow(M, 2)*pow(beta, 2) + 8*pow(M, 2)*beta + 2*pow(M, 2)) - 9*(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0)*(4*pow(M, 2)*pow(beta, 2)*pow(p0, 2) + 8*pow(M, 2)*beta*pow(p0, 2) + 4*pow(M, 2)*pow(p0, 2) - 16*M*beta*p0*q - 8*M*p0*q + 4*pow(beta, 2)*pow(p0, 2) + 4*beta*pow(p0, 2) + pow(p0, 2))/((4*M*pow(beta, 2) + 4*M*beta + M)*(16*pow(M, 2)*pow(beta, 2) + 16*pow(M, 2)*beta + 4*pow(M, 2))) + 2*pow(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0, 3)/pow(4*M*pow(beta, 2) + 4*M*beta + M, 3), 2)) + (27.0L/2.0L)*(-2*M*beta*pow(p0, 3) + 2*M*beta*pow(p0, 2)*p - 2*M*pow(p0, 2)*p + 2*M*p0*pow(p, 2) + 2*beta*pow(p0, 2)*q + pow(p0, 2)*q)/(8*pow(M, 2)*pow(beta, 2) + 8*pow(M, 2)*beta + 2*pow(M, 2)) - 9.0L/2.0L*(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0)*(4*pow(M, 2)*pow(beta, 2)*pow(p0, 2) + 8*pow(M, 2)*beta*pow(p0, 2) + 4*pow(M, 2)*pow(p0, 2) - 16*M*beta*p0*q - 8*M*p0*q + 4*pow(beta, 2)*pow(p0, 2) + 4*beta*pow(p0, 2) + pow(p0, 2))/((4*M*pow(beta, 2) + 4*M*beta + M)*(16*pow(M, 2)*pow(beta, 2) + 16*pow(M, 2)*beta + 4*pow(M, 2))) + pow(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0, 3)/pow(4*M*pow(beta, 2) + 4*M*beta + M, 3)) - 1.0L/3.0L*(-pow(M, 2)*pow(beta, 2)*p0 - 2*pow(M, 2)*beta*p0 - pow(M, 2)*p0 + 4*M*beta*q + 2*M*q - 4*pow(beta, 2)*p0 - 4*beta*p0 - p0)/(4*M*pow(beta, 2) + 4*M*beta + M);

        /////////////////////////////////////////////////////////////////////////

        T M_sq = M * M;
        T beta_sq = beta * beta;
        T p0_sq = p0 * p0;
        T p0_cu = p0_sq * p0;
        T p_sq = p * p;

        T tmp9 = (4 * M_sq * beta_sq * p0_sq + 8 * M_sq * beta * p0_sq + 4 * M_sq * p0_sq - 16 * M * beta * p0 * q - 8 * M * p0 * q + 4 * beta_sq * p0_sq + 4 * beta * p0_sq + p0_sq);
        T tmp10 = (16 * M_sq * beta_sq + 16 * M_sq * beta + 4 * M_sq);
        T tmp11 = -M_sq * beta_sq * p0 - 2 * M_sq * beta * p0 - M_sq * p0 + 4 * M * beta * q + 2 * M * q - 4 * beta_sq * p0 - 4 * beta * p0 - p0;
        T tmp12 = (-2 * M * beta * p0_cu + 2 * M * beta * p0_sq * p - 2 * M * p0_sq * p + 2 * M * p0 * p_sq + 2 * beta * p0_sq * q + p0_sq * q);
        T tmp13 = (8 * M_sq * beta_sq + 8 * M_sq * beta + 2 * M_sq);
        T tmp14 = (4 * M * beta_sq + 4 * M * beta + M);

        T tmp7 = tmp14 * tmp14;
        T tmp3 = tmp7 * tmp14;
        T tmp5 = tmp11 * tmp11;
        T tmp4 = tmp5 * tmp11;

        T tmp4_over_tmp3 = tmp4 / tmp3;
        T tmp9_over_tmp10 = tmp9 / tmp10;
        T tmp5_over_tmp7 = tmp5 / tmp7;
        T tmp12_over_tmp13 = tmp12 / tmp13;
        T tmp91410 = tmp9 / (tmp14 * tmp10);

        T tmp6_1 = 27 * tmp12_over_tmp13 - 9 * tmp11 * tmp91410 + 2 * tmp4_over_tmp3;
        T tmp6 = tmp6_1 * tmp6_1;
        T tmp8_1 = -3 * tmp9_over_tmp10 + tmp5_over_tmp7;
        T tmp8 = tmp8_1 * tmp8_1 * tmp8_1;

        T tmp2_inside = -4 * tmp8 + tmp6;
        if (tmp2_inside < 0) {
#ifdef WARNINGS
            debug("AnalQuadReturnMapping: WARNING Square root of negative number, tmp2_inside = ", tmp2_inside);
            debug("                   ... with these values: scale = ", scale);
            debug("                   ... with these values: p0    = ", p0, ", p0*scale = ", p0 * scale, " Pa");
            debug("                   ... with these values: p     = ", p, ", p*scale  = ", p * scale, " Pa");
            debug("                   ... with these values: q     = ", q, ", q*scale  = ", q * scale, " Pa");
            debug("          Compare to: ");
            debug("                  13.5*tmp12_over_tmp13 = ", 13.5 * tmp12_over_tmp13);
            debug("                  4.5*tmp11*tmp91410    = ", 4.5 * tmp11 * tmp91410);
            debug("                  tmp4_over_tmp3        = ", tmp4_over_tmp3);
#endif
            if (tmp2_inside < -1e-15) {
#ifdef WARNINGS
                debug("          tmp2_inside = ", tmp2_inside);
#endif
                T assumption_1 =                          13.5 * tmp12_over_tmp13 - 4.5 * tmp11 * tmp91410 + tmp4_over_tmp3;
                T assumption_2 = 0.5*sqrt(-tmp2_inside) + 13.5 * tmp12_over_tmp13 - 4.5 * tmp11 * tmp91410 + tmp4_over_tmp3;
                if ( abs(assumption_1-assumption_2)/abs(assumption_1) < 1e-2 ){
                    tmp2_inside = 0;
#ifdef WARNINGS
                    debug("                  The relative mistake by neglecting tmp2 is less than 1 percent");
#endif
                } else {
                    debug("                  FATAL The relative mistake by neglecting tmp2 is NOT less than 1 percent");
                    exit = 1;
                }
            } else {
                tmp2_inside = 0;
#ifdef WARNINGS
                debug("                  Even though the inside of the square root is negative, it is so small, on the order of machine precision");
#endif
            }
        }
        T tmp2 = sqrt(tmp2_inside);
        T tmp1 = cbrt(0.5 * tmp2 + 13.5 * tmp12_over_tmp13 - 4.5 * tmp11 * tmp91410 + tmp4_over_tmp3);

        if (tmp1 < 1e-15){
            debug("AnalQuadReturnMapping: FATAL tmp1 = ", tmp1);
            exit = 1;
        }

        T lamb = -0.3333333333333333333333333333333333333333333333333333333333333 * tmp8_1 / tmp1 - 0.3333333333333333333333333333333333333333333333333333333333333 * tmp1 - 0.3333333333333333333333333333333333333333333333333333333333333 * tmp11 / tmp14;

        /////////////////////////////////////////////////////////////////////////

        p = p0 * (-M * beta * lamb + M * lamb - p) / (2 * M * lamb - p0);
        q = beta * lamb + lamb / 2 + q;

        p0 *= scale;
        p *= scale;
        q *= scale;

        p = max(min(p, max_p), min_p);
        q = min(max(q, T(0)), max_q);

        T yield_function = q * (1 + 2 * beta) + 2 * M * (p + beta * p0) * (p - p0) / p0; // yield surface, do not use precomps here!
        if (std::abs(yield_function) > T(5)) {
#ifdef WARNINGS
            debug("AnalQuadReturnMapping: WARNING y(p,q) = ", yield_function, " at p = ", p, " and q = ", q);
#endif
            if (yield_function*yield_function > T(1e-4)*(p*p+q*q)) {
                debug("AnalQuadReturnMapping: FATAL y(p,q) = ", yield_function, " at p = ", p, " and q = ", q);
                exit = 1;
            }
        }
        if (!std::isfinite(p) || !std::isfinite(q)){
            debug("AnalQuadReturnMapping: FATAL p or q not finite");
            exit = 1;
        }
        return true;
    }
    return false;
}
