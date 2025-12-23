#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include "tool.h"
#include "gaussSieve.h"
#include <iostream>

using namespace std;
using namespace NTL;

int main(int argc, char *argv[])
{
    RR::SetPrecision(256);

    if (argc != 2)
    {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string input_file = argv[1];

    mat_ZZ B;
    if (!read_basis<ZZ>(input_file, B))
    {
        cerr << "Error: Failed to read basis from " << input_file << endl;
        return 1;
    }

    {
        ZZ det_unused;
        double delta = 0.99;
        long a = static_cast<long>(delta * 100);
        long b = 100;
        LLL(det_unused, B, a, b, 0);
    }

    cout << "Loaded basis:" << endl;
    cout << B << endl;
    RR gh = computeGH(B);
    cout << "Gaussian Heuristic: " << gh << endl;

    GaussSieveConfig config;
    config.sigma = to_RR("600");
    config.mu = 1.05;
    config.maxCycles = 1000;
    config.logEveryIteration = true;
    config.useGHStop = true;

    RR target_norm = to_RR(config.mu) * gh;
    cout << "Target mu: " << config.mu << endl;
    cout << "Target norm (mu * GH): " << target_norm << endl;

    vec_ZZ shortest = gaussSieveMain(B, config);

    cout << "Shortest vector (approx):" << endl;
    cout << shortest << endl;
    RR shortest_norm = Computenorm(shortest);
    cout << "Shortest norm: " << shortest_norm << endl;
    cout << "Shortest/GH ratio: " << (shortest_norm / gh) << endl;

    return 0;
}
