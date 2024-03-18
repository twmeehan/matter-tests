#ifndef PLASTIC_MODELS_HPP
#define PLASTIC_MODELS_HPP

#include "tools.hpp"

bool MCCRMA(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K, T rma_prefac);
bool MCCHardRMA(   T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T rma_prefac, T epv);
bool MCCHardExpRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T rma_prefac, T epv);
bool SinterMCCRMA(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T dt, T Sc, T tc, T ec, T epv, T S);

#endif // PLASTIC_MODELS_HPP
