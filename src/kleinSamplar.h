#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/BasicThreadPool.h>
#include <random>
#include <NTL/LLL.h>
#include "tool.h"

using namespace NTL;
using namespace std;

ZZ sample_gaussian_1D(const RR &center, const RR &sigma);
vec_ZZ KleinSampler(const mat_ZZ &B, const RR &sigma, const vec_RR &c_target);