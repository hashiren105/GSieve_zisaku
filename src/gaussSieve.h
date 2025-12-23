#ifndef GAUSS_SIEVE_H
#define GAUSS_SIEVE_H

#include <vector>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/RR.h>
#include "tool.h"

NTL_CLIENT

struct GaussSieveConfig {
    RR sigma;
    double mu;
    int maxCycles = 100;
    bool logEveryIteration = true;
    bool useGHStop = false;
};

vec_ZZ gaussReduce(const vec_ZZ &p, std::vector<vec_ZZ> &L, std::vector<vec_ZZ> &S);
vec_ZZ gaussSieveMain(const mat_ZZ &B, const GaussSieveConfig &config);

#endif // GAUSS_SIEVE_H
