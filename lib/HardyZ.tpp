#include "Theta.h"  
#include "Bernoulli.h"  
#include <cmath>
#include <vector>
#include <complex>
#include <numbers>      
#include <ranges>       
#include <algorithm>    

namespace Zeta::Hardy {

    namespace detail {

        template <std::floating_point T>
        std::complex<T> zetaEM(std::complex<T> s, int N) {
            if (N <= 1) return { T{0}, T{0} };

            auto range = std::views::iota(1, N);
            
            std::complex<T> sum = std::ranges::fold_left(
                range,
                std::complex<T>{0, 0},
                [s](std::complex<T> acc, int n) {
                    return acc + std::pow(static_cast<T>(n), -s);
                }
            );

            const T N_dbl = static_cast<T>(N);
            const std::complex<T> N_pow_minus_s = std::pow(N_dbl, -s);
            const T inv_N = T{1} / N_dbl;
            const T inv_N_sq = inv_N * inv_N;

            std::complex<T> term_integral = (N_dbl * N_pow_minus_s) / (s - T{1});
            std::complex<T> term_half     = T{0.5} * N_pow_minus_s;

            T Bern2 = Zeta::bernoulli<T>(2);
            T Bern4 = Zeta::bernoulli<T>(4);

            std::complex<T> term_B2 = (Bern2 / T{2}) * (-s) * (N_pow_minus_s * inv_N);
            
            std::complex<T> term_B4 = (Bern4 / T{24}) 
                                    * (-s) * (-s - T{1}) * (-s - T{2}) 
                                    * (N_pow_minus_s * inv_N * inv_N_sq);

            return sum + term_integral + term_half + term_B2 + term_B4;
        }

        template <std::floating_point T>
        T computeEM(T t) {
            const int N = std::max(static_cast<int>(std::abs(t)) + 5, 15);

            std::complex<T> s(T{0.5}, t);
            std::complex<T> zeta_val = zetaEM(s, N);
            
            T theta_val = Zeta::theta<T>(t); 

            std::complex<T> phase(T{0}, theta_val);
            return (std::exp(phase) * zeta_val).real();
        }

        template <std::floating_point T>
        T computeRS(T t) {
            constexpr T PI = std::numbers::pi_v<T>;
            
            int N = static_cast<int>(std::floor(std::sqrt(t / (T{2} * PI))));
            if (N < 1) return T{0};

            T theta_val = Zeta::theta<T>(t);

            // Formula: Sum[ cos(theta - t*ln(n)) / sqrt(n) ]
            auto range = std::views::iota(1, N + 1);

            T sum = std::ranges::fold_left(
                range,
                T{0},
                [theta_val, t](T acc, int n) {
                    T n_val = static_cast<T>(n);
                    T term_arg = theta_val - (t * std::log(n_val));
                    return acc + (std::cos(term_arg) / std::sqrt(n_val));
                }
            );

            return T{2} * sum;
        }

        template <std::floating_point T>
        std::vector<T> computeOS(T start_t, T length, int points) {
            constexpr T PI = std::numbers::pi_v<T>;

            std::vector<T> results;
            results.reserve(points);

            int N = static_cast<int>(std::floor(std::sqrt(start_t / (T{2} * PI))));
            if (N < 1) N = 1;

            std::vector<std::complex<T>> base_terms(N + 1);
            
            std::ranges::for_each(
                std::views::iota(1, N + 1),
                [start_t, &base_terms](int n) {
                    T n_val = static_cast<T>(n);
                    T ln_n = std::log(n_val);
                    T mag = T{1} / std::sqrt(n_val);
                    T phase = -start_t * ln_n;
                    
                    base_terms[n] = std::polar(mag, phase);
                }
            );

            T step = (points > 1) ? (length / static_cast<T>(points - 1)) : T{0};
            auto indices = std::views::iota(1, N + 1); 

            return std::views::iota(0, points) 
                    | std::views::transform(
                    [start_t, step, indices, &base_terms](int k) {
                        T delta = static_cast<T>(k) * step;
                        T t_current = start_t + delta;

                        T theta_val = Zeta::theta<T>(t_current);
                        std::complex<T> rot_phase = std::polar(T{1}, theta_val);

                        std::complex<T> sum = std::ranges::fold_left(
                            indices,
                            std::complex<T>{0, 0},
                            [delta, &base_terms](std::complex<T> acc, int n) {
                                T ln_n = std::log(static_cast<T>(n));
                                T perturbation_phase = -delta * ln_n;
                                
                                std::complex<T> term = base_terms[n] * std::polar(T{1}, perturbation_phase);
                                return acc + term;
                            }
                        );

                        return T{2} * (rot_phase * sum).real();
                    }
                ) 
                | std::ranges::to<std::vector<T>>(); 
        }

    } 

    // =====================================================================
    // Public API Implementation
    // =====================================================================

    template <std::floating_point T>
    T compute(T t, Method method) {
        if (std::abs(t) < T{1e-9}) return T{-0.5}; 

        switch (method) {
            case Method::EulerMaclaurin:
                return detail::computeEM<T>(t);
            case Method::RiemannSiegel:
                return detail::computeRS<T>(t);
            default:
                return T{0};
        }
    }

    template <std::floating_point T>
    std::vector<T> computeBlock(T start_t, T length, int points, Method method) {
        if (method == Method::OdlyzkoSchonhage) {
            return detail::computeOS<T>(start_t, length, points);
        }

        T step = (points > 1) ? (length / static_cast<T>(points - 1)) : T{0};
        auto indices = std::views::iota(0, points);

        return indices 
            | std::views::transform([start_t, step, method](int i) {
                T t = start_t + (static_cast<T>(i) * step);
                return compute<T>(t, method);
            })
            | std::ranges::to<std::vector<T>>();
    }

}