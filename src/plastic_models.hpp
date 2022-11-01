#ifndef PLASTIC_MODELS_HPP
#define PLASTIC_MODELS_HPP

#include "tools.hpp"

bool MCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K);
bool MCCHardRMA(T& p, T& q, int& exit, T M, T beta, T mu, T K, T xi, T epv);
bool MCCHardExpRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T epv);
bool PerzynaMCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T dt, T d, T perzyna_visc);
bool PerzynaMCCHardRMA(T& p, T& q, int& exit, T M, T p00, T beta, T xi, T mu, T K, T dt, T d, T perzyna_visc, T epv);
bool SinterMCCRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T dt, T Sc, T tc, T ec, T epv, T S);

#endif // PLASTIC_MODELS_HPP
