#ifndef KLEIN_SAMPLAR_H
#define KLEIN_SAMPLAR_H

#include <NTL/mat_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/vec_ZZ.h>

using namespace NTL;

ZZ sample_gaussian_1D(const RR &center, const RR &sigma);
vec_ZZ KleinSampler(const mat_ZZ &B, const RR &sigma, const vec_RR &c_target);
bool KleinSamplerAborted();

#endif // KLEIN_SAMPLAR_H
