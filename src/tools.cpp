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

bool CamClayReturnMapping(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus)
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

bool QuadraticReturnMapping(T& p, T& q, int& exit, T trace_epsilon, T norm_eps_hat, T M, T p0, T beta, T mu, T bulk_modulus)
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
        for (int iter = 0; iter < 400; iter++) {

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
        }
        p0 *= scale;
        p *= scale;
        q *= scale;

        p = max(min(p, max_p), min_p);
        q = min(abs(q), max_q);
        assert((q * (1 + 2 * beta) + 2 * M * (p + beta * p0) * (p - p0) / p0) <= T(1e-3)); // yield surface, do not use precomps here!
        assert(std::isfinite(p));
        assert(std::isfinite(q));
        return true;
    }
    return false;
}



bool AnalQuadReturnMapping(T& p, T& q, int& exit, T M, T p0, T beta)
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
            debug("AnalQuadReturnMapping: Square root of negative number, tmp2_inside = ", tmp2_inside);
            debug("                   ... with these values: scale = ", scale);
            debug("                   ... with these values: p0    = ", p0, ", p0*scale = ", p0 * scale, " Pa");
            debug("                   ... with these values: p     = ", p, ", p*scale  = ", p * scale, " Pa");
            debug("                   ... with these values: q     = ", q, ", q*scale  = ", q * scale, " Pa");
            debug("          Compare to: ");
            debug("                  13.5*tmp12_over_tmp13 = ", 13.5 * tmp12_over_tmp13);
            debug("                  4.5*tmp11*tmp91410    = ", 4.5 * tmp11 * tmp91410);
            debug("                  tmp4_over_tmp3        = ", tmp4_over_tmp3);
            if (tmp2_inside < -1e-15) {
                debug("          tmp2_inside = ", tmp2_inside);
                T assumption_1 =                          13.5 * tmp12_over_tmp13 - 4.5 * tmp11 * tmp91410 + tmp4_over_tmp3;
                T assumption_2 = 0.5*sqrt(-tmp2_inside) + 13.5 * tmp12_over_tmp13 - 4.5 * tmp11 * tmp91410 + tmp4_over_tmp3;
                if ( abs(assumption_1-assumption_2)/abs(assumption_1) < 1e-2 ){
                    tmp2_inside = 0;
                    debug("                  The relative mistake by neglecting tmp2 is less than 1 percent");
                } else {
                    debug("                  The relative mistake by neglecting tmp2 is NOT less than 1 percent");
                    exit = 1;
                }
            } else {
                tmp2_inside = 0;
                debug("                  Even though the inside of the square root is negative, it is so small, on the order of machine precision");
            }
        }
        T tmp2 = sqrt(tmp2_inside);
        T tmp1 = cbrt(0.5 * tmp2 + 13.5 * tmp12_over_tmp13 - 4.5 * tmp11 * tmp91410 + tmp4_over_tmp3);

        if (tmp1 < 1e-15){
            debug("AnalQuadReturnMapping: tmp1 = ", tmp1);
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
        if (yield_function > T(1)) {
            debug("AnalQuadReturnMapping: yield_function = ", yield_function);
            exit = 1;
        }
        if (!std::isfinite(p) || !std::isfinite(q)){
            debug("AnalQuadReturnMapping: p or q not finite");
            exit = 1;
        }
        return true;
    }
    return false;
}
