#include "kleinSamplar.h"
#include "tool.h"
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_RR.h>
#include <NTL/vec_ZZ.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

bool g_sampler_aborted = false;

RR random_unit_interval()
{
    static const long UNIFORM_BITS = 256;
    static const RR denom = [] {
        ZZ one = ZZ(1);
        one <<= UNIFORM_BITS;
        return to_RR(one);
    }();

    while (true) {
        ZZ sample;
        RandomLen(sample, UNIFORM_BITS);
        if (IsZero(sample))
            continue;
        RR value = to_RR(sample) / denom;
        if (value > 0 && value < 1)
            return value;
    }
}

RR standard_gaussian()
{
    static const RR two_pi = 2 * ComputePi_RR();

    while (true) {
        RR u1 = random_unit_interval();
        if (IsZero(u1))
            continue;
        RR u2 = random_unit_interval();
        RR radius = sqrt(-2 * log(u1));
        RR angle = two_pi * u2;
        RR z;
        cos(z, angle);
        z *= radius;
        return z;
    }
}

mat_RR computeGramSchmidt(const mat_ZZ& B)
{
    long n = B.NumRows();
    long dim = B.NumCols();

    mat_RR G;
    G.SetDims(n, dim);

    for (long i = 0; i < n; ++i)
    {
        vec_RR bi = to_RR(B[i]);

        for (long j = 0; j < i; ++j)
        {
            vec_RR bj_star = G[j];
            RR denom;
            InnerProduct(denom, bj_star, bj_star);
            if (IsZero(denom) || denom <= 0)
                continue;

            RR mu;
            InnerProduct(mu, bi, bj_star);
            mu /= denom;
            bi -= mu * bj_star;
        }

        G[i] = bi;
    }

    return G;
}

std::string timestamp()
{
    using namespace std::chrono;
    auto now = system_clock::now();
    std::time_t now_t = system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&now_t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%F %T");
    return oss.str();
}

void log_message(const std::string &message)
{
    std::cout << "[" << timestamp() << "] [KleinSampler] " << message << std::endl;
}

} // namespace

bool KleinSamplerAborted()
{
    return g_sampler_aborted;
}

ZZ sample_gaussian_1D(const RR &center, const RR &sigma)
{
    if (!(sigma > 0))
        return RoundToZZ(center);

    RR sample = center + sigma * standard_gaussian();
    return RoundToZZ(sample);
}

vec_ZZ KleinSampler(const mat_ZZ &B, const RR &sigma, const vec_RR &c_target)
{
    g_sampler_aborted = false;
    int n = B.NumRows();
    int dim = B.NumCols();

    mat_RR B_star = computeGramSchmidt(B);

    const long maxAttempts = 1000;
    long attempt = 0;
    for (;;)
    {
        ++attempt;
        vec_ZZ v;
        v.SetLength(dim);
        clear(v);

        vec_RR c_vec = c_target;

        for (int i = n - 1; i >= 0; --i)
        {
            vec_ZZ b_i = B[i];
            vec_RR b_star_i = B_star[i];

            RR b_i_norm_sq;
            InnerProduct(b_i_norm_sq, b_star_i, b_star_i);
            if (IsZero(b_i_norm_sq) || b_i_norm_sq <= 0)
                continue;

            RR dot;
            InnerProduct(dot, c_vec, b_star_i);
            RR d_i = dot / b_i_norm_sq;

            RR sigma_i = sigma / sqrt(b_i_norm_sq);
            ZZ z_i = sample_gaussian_1D(d_i, sigma_i);

            vec_RR b_i_RR = to_RR(b_i);
            c_vec -= to_RR(z_i) * b_i_RR;
            v += z_i * b_i;
        }

        if (!IsZero(v)) {
            std::ostringstream status;
            status << "Sampled non-zero vector after " << attempt
                   << " attempt(s) | norm=" << Computenorm(v);
            log_message(status.str());
            return v;
        }
        std::ostringstream status;
        status << "Sampling failed (zero vector), attempt=" << attempt;
        log_message(status.str());
        if (attempt >= maxAttempts) {
            log_message("Sampling aborted after maxAttempts=" + std::to_string(maxAttempts));
            g_sampler_aborted = true;
            vec_ZZ zero;
            zero.SetLength(dim);
            clear(zero);
            return zero;
        }
    }
}
