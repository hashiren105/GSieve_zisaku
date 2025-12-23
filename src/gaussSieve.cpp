#include "gaussSieve.h"
#include "kleinSamplar.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

RR norm_squared_rr(const vec_ZZ &v)
{
    vec_RR r = to_RR(v);
    RR result;
    InnerProduct(result, r, r);
    return result;
}

vec_ZZ subtract_vectors(const vec_ZZ &a, const vec_ZZ &b)
{
    long len = std::max(a.length(), b.length());
    vec_ZZ result;
    result.SetLength(len);
    for (long i = 0; i < len; ++i) {
        ZZ ai = (i < a.length()) ? a[i] : ZZ(0);
        ZZ bi = (i < b.length()) ? b[i] : ZZ(0);
        result[i] = ai - bi;
    }
    return result;
}

bool is_zero_vector(const vec_ZZ &v)
{
    for (long i = 0; i < v.length(); ++i) {
        if (!IsZero(v[i]))
            return false;
    }
    return true;
}

std::string rr_to_string(const RR &value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
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
    std::cout << "[" << timestamp() << "] [GaussSieve] " << message << std::endl;
}

} // namespace

vec_ZZ gaussReduce(const vec_ZZ &p_in, std::vector<vec_ZZ> &L, std::vector<vec_ZZ> &S)
{
    vec_ZZ p = p_in;
    RR norm_p = norm_squared_rr(p);
    bool updated = true;

    while (updated) {
        updated = false;
        for (size_t i = 0; i < L.size(); ++i) {
            vec_ZZ vi = L[i];
            RR norm_vi = norm_squared_rr(vi);
            vec_ZZ p_minus_vi = subtract_vectors(p, vi);
            RR norm_diff = norm_squared_rr(p_minus_vi);

            if (norm_vi < norm_p && norm_diff < norm_p) {
                p = p_minus_vi;
                norm_p = norm_diff;
                updated = true;
                break;
            }
        }
    }

    if (is_zero_vector(p)) {
        return p;
    }

    for (size_t i = 0; i < L.size();) {
        vec_ZZ vi = L[i];
        RR norm_vi = norm_squared_rr(vi);
        vec_ZZ vi_minus_p = subtract_vectors(vi, p);
        RR norm_diff = norm_squared_rr(vi_minus_p);

        if (norm_vi > norm_p && norm_diff < norm_vi) {
            if (i != L.size() - 1) {
                L[i] = L.back();
            }
            L.pop_back();
            if (!is_zero_vector(vi_minus_p)) {
                S.push_back(vi_minus_p);
            }
        } else {
            ++i;
        }
    }

    return p;
}

vec_ZZ gaussSieveMain(const mat_ZZ &B, const GaussSieveConfig &config)
{
    std::vector<vec_ZZ> L;
    std::vector<vec_ZZ> S;
    int c = config.maxCycles;
    int K = 0;
    RR gh;
    RR threshold;
    if (config.useGHStop) {
        gh = computeGH(B);
        threshold = to_RR(config.mu) * gh;
    }
    vec_RR c_target;
    c_target.SetLength(B.NumCols());
    clear(c_target);

    vec_ZZ v_new;
    int iteration = 0;
    long long collisionCount = 0;

    {
        std::ostringstream oss;
        oss << "Starting sweep (sigma=" << rr_to_string(config.sigma) << ", mu=" << config.mu
            << ", maxCycles=" << c << ", GHstop=" << (config.useGHStop ? "on" : "off") << ")";
        log_message(oss.str());
    }

    while (K < c) {
        ++iteration;
        if (!S.empty()) {
            v_new = S.back();
            S.pop_back();
        } else {
            if (config.logEveryIteration) {
                log_message("Sampling from KleinSampler");
            }
            v_new = KleinSampler(B, config.sigma, c_target);
            if (KleinSamplerAborted()) {
                log_message("Aborting GaussSieve: KleinSampler hit maxAttempts");
                break;
            }
        }

        if (is_zero_vector(v_new)) {
            ++K;
            if (config.logEveryIteration) {
                std::ostringstream status;
                status << "Zero vector encountered; K=" << K << " / " << c;
                log_message(status.str());
            }
            continue;
        }

        v_new = gaussReduce(v_new, L, S);

        if (is_zero_vector(v_new)) {
            ++K;
            ++collisionCount;
            if (config.logEveryIteration) {
                std::ostringstream status;
                status << "Collision (reduced to zero); collisions=" << collisionCount
                       << ", K=" << K << " / " << c;
                log_message(status.str());
            }
            continue;
        }

        RR norm_v = Computenorm(v_new);

        if (config.logEveryIteration) {
            std::ostringstream status;
            status << "Iteration " << iteration << " | norm=" << rr_to_string(norm_v)
                   << " | L=" << L.size() << " | S=" << S.size();
            log_message(status.str());
        }

        if (config.useGHStop && norm_v <= threshold) {
            log_message("Terminating early: norm <= " + std::to_string(config.mu) + " × GH");
            break;
        }

        L.push_back(v_new);
        if (config.logEveryIteration) {
            std::ostringstream status;
            status << "Added vector to L, size=" << L.size();
            log_message(status.str());
        }
    }

    return v_new;
}
