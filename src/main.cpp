#include <NTL/mat_ZZ.h>
#include "tool.h"
#include "kleinSamplar.h"
#include <iostream>

using namespace std;
using namespace NTL;

int main(int argc, char *argv[])
{
    RR::SetPrecision(256);

    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    string input_file = argv[1];
    string output_file = argv[2];

    mat_ZZ B;

    // Read basis using tool.cpp function
    if (!read_basis<ZZ>(input_file, B))
    {
        cerr << "Error: Failed to read basis from " << input_file << endl;
        return 1;
    }

    // Output basis to stdout
    cout << "Loaded basis:" << endl;
    cout << B << endl;

    RR gh = computeGH(B);
    cout << "Gaussian Heuristic: " << gh << endl;

    // σ ≈ 0.3325765549e121 (n = 40)
    RR sigma = to_RR("0.3325765549e121");
    vec_RR c_target;
    c_target.SetLength(B.NumCols());
    for (long i = 0; i < c_target.length(); ++i)
    {
        c_target[i] = 0;
    }

    vec_ZZ v = KleinSampler(B, sigma, c_target);

    cout << "Sampled vector v:" << endl;
    cout << v << endl;

    return 0;
}

