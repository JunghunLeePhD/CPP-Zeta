#include "HardyZ.h"
#include "Theta.h"  
#include "Bernoulli.h"
#include <cmath>
#include <vector>
#include <iostream>
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

        double Bern2 = Zeta::Bernoulli::get(2);
        double Bern4 = Zeta::Bernoulli::get(4);

        std::complex<double> term_B2 = (Bern2 / 2.0) * (-s) * (N_pow_minus_s * inv_N);
        std::complex<double> term_B4 = (Bern4 / 24.0) * (-s) * (-s - 1.0) * (-s - 2.0) * (N_pow_minus_s * inv_N * inv_N_sq);

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

    double HardyZ::computeRS(double t) {
        int N = static_cast<int>(std::floor(std::sqrt(t / (2.0 * M_PI))));
        if (N <= 1) return 0.0;

        double theta = Theta::value(t);

        std::vector<int> indices(N);
        std::iota(indices.begin(), indices.end(), 1);

        double sum = std::accumulate(
            indices.begin(), indices.end(),
            0.0,
            [theta, t](double current_sum, int n) {
                double n_dbl = static_cast<double>(n);
                double term_arg = theta - (t * std::log(n_dbl));
                return current_sum + (std::cos(term_arg) / std::sqrt(n_dbl));
            }
        );

        return 2.0 * sum;
    }

    std::vector<double> HardyZ::computeOS(double start_t, double length, int points) {
        std::vector<double> results;
        results.reserve(points);

        int N = static_cast<int>(std::floor(std::sqrt(start_t / (2.0 * M_PI))));
        if (N < 1) N = 1;

        std::vector<std::complex<double>> base_terms(N + 1);
        
        for (int n = 1; n <= N; ++n) {
            double ln_n = std::log(static_cast<double>(n));
            double mag = 1.0 / std::sqrt(static_cast<double>(n));
            double phase = -start_t * ln_n;
            
            base_terms[n] = std::polar(mag, phase);
        }
        std::vector<int> indices(N);
        std::iota(indices.begin(), indices.end(), 1);

        double step = (points > 1) ? (length / (points - 1)) : 0.0;

        for (int k = 0; k < points; ++k) {
            double delta = k * step;       
            double t_current = start_t + delta;

            double theta = Theta::value(t_current);
            std::complex<double> rot_phase = std::polar(1.0, theta);

            std::complex<double> sum = std::accumulate(
                indices.begin(), indices.end(), 
                std::complex<double>(0.0),
                [delta, &base_terms](std::complex<double> current_sum, int n) { 
                    double ln_n = std::log(static_cast<double>(n));
                    double perturbation_phase = -delta * ln_n;
                    
                    std::complex<double> term = base_terms[n] * std::polar(1.0, perturbation_phase);

                    return current_sum + term; 
                }
            );

            double Z_val = 2.0 * (rot_phase * sum).real();
            
            results.push_back(Z_val);
        }

        return results;
    }

    double HardyZ::compute(double t, Method method) {
        if (std::abs(t) < 1e-9) return -0.5;

        switch (method) {
            case Method::EulerMaclaurin:
                return computeEM(t);
            case Method::RiemannSiegel:
                return computeRS(t);
            default:
                return 0;
        }
    }

    std::vector<double> HardyZ::computeBlock(double start_t, double length, int points, Method method) {
        if (method == Method::OdlyzkoSchonhage) {
            return computeOS(start_t, length, points);
        }

        std::vector<double> results;
        results.reserve(points);

        double step = (points > 1) ? (length / (points - 1)) : 0.0;
        
        for (int i = 0; i < points; ++i) {
            double t = start_t + (i * step);
            results.push_back(compute(t, method));
        }

        return results;
    }
}