#include "HardyZ.h"
#include "Theta.h"  
#include <cmath>
#include <vector>
#include <iostream>
#include <numeric>
#include <execution> 

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Zeta {

    std::complex<double> HardyZ::zetaEM(std::complex<double> s, int N) {
        if (N <= 1) return 0.0;

        std::vector<int> indices(N - 1);
        std::iota(indices.begin(), indices.end(), 1);

        std::complex<double> sum = std::accumulate(
            indices.begin(), indices.end(),
            std::complex<double>(0.0), 
            [s](std::complex<double> current_sum, int n) {
                return current_sum + std::pow(static_cast<double>(n), -s);
            }
        );

        const double N_dbl = static_cast<double>(N);
        
        const std::complex<double> N_pow_minus_s = std::pow(N_dbl, -s);
        const double inv_N = 1.0 / N_dbl;
        double inv_N_sq = inv_N * inv_N;
        
        std::complex<double> term_integral = (N_dbl * N_pow_minus_s) / (s - 1.0);
        std::complex<double> term_half = 0.5 * N_pow_minus_s;

        std::complex<double> term_B2 = (1.0 / 12.0) * s * (N_pow_minus_s * inv_N);
        std::complex<double> term_B4 = (-1.0 / 720.0) * s * (s + 1.0) * (s + 2.0) * (N_pow_minus_s * inv_N * inv_N_sq);

        return sum + term_integral + term_half + term_B2 + term_B4;
    }

    double HardyZ::computeEM(double t) {
        int N = static_cast<int>(std::abs(t)) + 5;
        if (N < 15) N = 15;

        std::complex<double> s(0.5, t);
        std::complex<double> zeta = zetaEM(s, N);
        
        double theta = Theta::value(t); 

        std::complex<double> phase(0.0, theta);
        return (std::exp(phase) * zeta).real();
    }

    double HardyZ::compute(double t, Method method) {
        if (std::abs(t) < 1e-9) return -0.5;

        switch (method) {
            case Method::EulerMaclaurin:
                return computeEM(t);
            default:
                return 0;
        }
    }
}