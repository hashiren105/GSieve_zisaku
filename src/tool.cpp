#include "tool.h"

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

    ZZ detB_ZZ = determinant(B);
    RR detB = to_RR(detB_ZZ);
    RR pi = ComputePi_RR();
    RR e = exp(to_RR(1));

    RR coeff = sqrt(to_RR(n) / (2 * pi * e));
    RR vol_root = pow(abs(detB), to_RR(1.0) / n);

    return coeff * vol_root;
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