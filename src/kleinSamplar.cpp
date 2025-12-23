#include "kleinSamplar.h"

ZZ sample_gaussian_1D(const RR &center, const RR &sigma)
{
    static random_device rd;
    static mt19937 gen(rd());

    double c = to_double(center);
    double s = to_double(sigma);

    normal_distribution<double> dist(c, s);
    double sampled = dist(gen);

    RR rr_sampled = to_RR(sampled);
    RR rr_rounded = round(rr_sampled);

    return RoundToZZ(rr_rounded);
}

vec_ZZ KleinSampler(const mat_ZZ &B, const RR &sigma, const vec_RR &c_target)
{
    int n = B.NumRows();
    
    int dim = B.NumCols();

    // 初期化
    vector<vec_ZZ> v_list;
    v_list.resize(n);
    for (int i = 0; i < n; ++i)
    {
        v_list[i].SetLength(dim);
        clear(v_list[i]);
    }
    vec_RR c_vec;
    vector<vec_RR> c_list;
    c_list.resize(n);
    for (int i = 0; i < n; ++i)
    {
        c_list[i].SetLength(dim);
        clear(c_list[i]);
    }
    c_list[n - 1] = c_target;
    mat_ZZ B_star_ZZ = B;

    // GSO用データ構造
    mat_RR mu;
    vec_RR c;
    ComputeGS(B_star_ZZ, mu, c);
    mat_RR B_star = to_RR(B_star_ZZ);

    // メインループ
    for (int i = n - 1; i > 0; --i)
    {
        vec_ZZ b_i = B[i];
        vec_RR b_star_i = B_star[i];
        vec_RR c_vec = c_list[i];

        RR b_i_norm_sq;
        InnerProduct(b_i_norm_sq, b_star_i, b_star_i);

        RR dot;
        InnerProduct(dot, c_vec, b_star_i); // (c_i, b*_i)

        RR d_i = dot / b_i_norm_sq;

        RR sigma_i = sigma / sqrt(b_i_norm_sq);

        ZZ z_i = sample_gaussian_1D(d_i, sigma_i);

        vec_RR b_i_RR = to_RR(b_i);

        c_list[i-1] = c_list[i] - to_RR(z_i) * b_i_RR;
        v_list[i-1] = v_list[i] + z_i * b_i;
    }

    return v_list[0]; // Or adjust this return depending on final desired v
}

