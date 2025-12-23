#include "tool.h"
#include <cmath>

vec_RR to_RR(const vec_ZZ &v)
{
    vec_RR r;
    r.SetLength(v.length());
    for (long i = 0; i < v.length(); ++i)
        r[i] = to_RR(v[i]);
    return r;
}

mat_RR to_RR(const mat_ZZ &M)
{
    mat_RR R;
    R.SetDims(M.NumRows(), M.NumCols());
    for (long i = 0; i < M.NumRows(); ++i)
        for (long j = 0; j < M.NumCols(); ++j)
            R[i][j] = to_RR(M[i][j]);
    return R;
}

RR computeGH(const mat_ZZ &B)
{
    long n = B.NumRows();
    if (n == 0) {
        return RR(0);
    }

    mat_ZZ Bt;
    transpose(Bt, B);
    mat_ZZ G;
    mul(G, B, Bt);
    ZZ detG_ZZ = determinant(G);
    if (IsZero(detG_ZZ)) {
        return RR(0);
    }

    RR detG = abs(to_RR(detG_ZZ));
    RR vol = sqrt(detG);

    double nd = static_cast<double>(n);
    double pi = std::acos(-1.0);
    double logVn = (nd / 2.0) * std::log(pi) - std::lgamma(nd / 2.0 + 1.0);

    RR log_vol = log(vol);
    RR log_gh = (log_vol - to_RR(logVn)) / to_RR(n);
    return exp(log_gh);
}

bool Reduce(vec_ZZ *p1, const vec_ZZ &p2)
{
    if (Computenorm(*p1) < Computenorm(p2))
    {
        ZZ dot1, dot2;
        InnerProduct(dot1, *p1, p2);
        InnerProduct(dot2, p2, p2);

        RR a = to_RR(dot1);
        a *= 2;
        RR b = to_RR(dot2);

        RR abs_a = abs(a);
        if (abs_a > b)
        {
            RR q = a / b;
            ZZ q_round = RoundToZZ(q);
            *p1 -= q_round * p2;
            return true;
        }
    }
    return false;
}

RR Computenorm(const vec_ZZ &v)
{
    vec_RR v_RR = to_RR(v);
    RR norm;
    InnerProduct(norm, v_RR, v_RR);
    return sqrt(norm);
}
