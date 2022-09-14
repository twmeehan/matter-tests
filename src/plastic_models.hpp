#ifndef PLASTIC_MODELS_HPP
#define PLASTIC_MODELS_HPP

#include "tools.hpp"

bool ModifiedCamClayRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K);
bool ModifiedCamClayHardRMA(T& p, T& q, int& exit, T M, T epv, T beta, T mu, T K, T xi);
bool PerzynaMCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc);
bool PerzynaMCCHardRMA(T& p, T& q, int& exit, T M, T p00, T beta, T xi, T mu, T K, T dt, T d, T perzyna_visc, T epv);
bool PerzynaSinterMCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T epv, T S, T visc, T Sinf, T tc, T ec, T xi);

#endif // PLASTIC_MODELS_HPP
