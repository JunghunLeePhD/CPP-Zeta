#pragma once

#include <complex>
#include <vector>
#include <concepts>

namespace Zeta {

    /**
     * @brief Methods available for computing the Hardy Z function.
     */
    enum class Method {
        /**
         * @brief Euler-Maclaurin Summation.
         * * **Complexity:** $ O(t) $
         * * **Precision:** High (uses Bernoulli correction terms).
         * * **Use Case:** Recommended for $ t < 10,000 $ or high accuracy checks.
         */
        EulerMaclaurin,

        /**
         * @brief Riemann-Siegel Formula (Main Sum approximation).
         * * **Complexity:** $ O(\sqrt{t}) $
         * * **Precision:** Moderate (Main sum only, ignores $ \Psi $ remainders).
         * * **Use Case:** Recommended for large $ t $ (e.g., $ t > 10^5 $).
         */
        RiemannSiegel,

        /**
         * @brief Odlyzko-Schönhage Algorithm.
         * * **Complexity:** $ O(t^{1/3}) $ (amortized over a block).
         * * **Precision:** High.
         * * **Use Case:** Only efficient when computing **blocks** of zeros for extremely large $ t $.
         */
        OdlyzkoSchonhage
    };

    /**
     * @namespace Hardy
     * @brief Functions for computing the Hardy Z-function on the critical line.
     * * Definition: $ Z(t) = e^{i \theta(t)} \zeta\left(\frac{1}{2} + it\right) $
     */
    namespace Hardy {

        // =============================================================
        // Public API
        // =============================================================

        /**
         * @brief Computes the value of Z(t).
         * @tparam T Floating point type (float, double, long double).
         * @param t The imaginary component of the argument.
         * @param method The algorithm to use (default: EulerMaclaurin).
         * @return The real value $ Z(t) $.
         */
        template <std::floating_point T>
        [[nodiscard]]
        T compute(T t, Method method = Method::EulerMaclaurin);

        /**
         * @brief Computes a range of Z values efficiently.
         * * Primary entry point for Odlyzko-Schönhage blocks.
         * @param start_t The starting height.
         * @param length The length of the interval.
         * @param points The number of sampling points.
         * @param method Algorithm (default: OdlyzkoSchonhage).
         * @return A vector containing the Z values.
         */
        template <std::floating_point T>
        [[nodiscard]]
        std::vector<T> computeBlock(T start_t, T length, int points, Method method = Method::OdlyzkoSchonhage);

        namespace detail {
            
        /**
             * @brief Helper for computeEM to calculate the complex Zeta value.
             * Computes $ \zeta(s) $ via:
             * $$
             * \zeta(s) \approx \sum_{n=1}^{N-1} n^{-s} + \frac{N^{1-s}}{s-1} + \frac{1}{2}N^{-s} + \sum_{k=1}^{m} \frac{B_{2k}}{(2k)!} f^{(2k-1)}(N)
             * $$
             * @param s The complex argument $ \frac{1}{2} + it $.
             * @param N The summation cutoff limit.
             */
            template <std::floating_point T>
            std::complex<T> zetaEM(std::complex<T> s, int N);


            /**
             * @brief Computes Z(t) using the Euler-Maclaurin summation for $ \zeta(s) $.
             * Computes $ \zeta(s) $ via zetaEM 
             * and rotates the result by $ e^{i\theta(t)} $.
             */
            template <std::floating_point T>
            T computeEM(T t);

            /**
             * @brief Computes Z(t) using the Riemann-Siegel Main Sum.
             * Note that it does not use $ \zeta(s) $ as the Euler-Maclaurin method did.
             * Uses the approximation:
             * $$
             * Z(t) \approx 2 \sum_{n=1}^{\lfloor \sqrt{t/2\pi} \rfloor} \frac{\cos(\theta(t) - t \ln n)}{\sqrt{n}}
             * $$
             */
            template <std::floating_point T>
            T computeRS(T t);

            /**
             * @brief Implementation of the rational function expansion / FFT method.
             * * Uses the Taylor expansion of $ n^{-i\delta} $ to evaluate sums quickly.
             */
            template <std::floating_point T>
            std::vector<T> computeOS(T start_t, T length, int points);

        } 

    }
}

#include "HardyZ.tpp"