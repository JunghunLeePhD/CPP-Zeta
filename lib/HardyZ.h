#ifndef HARDYZ_H
#define HARDYZ_H

#include <complex>
#include <vector>

namespace Zeta {

    /**
     * @brief Methods available for computing the Hardy Z function.
     */
    enum class Method {
        /**
         * @brief Euler-Maclaurin Summation.
         * * **Complexity:** $ O(t) $
         * * **Precision:** High (uses Bernoulli correction terms).
         * * **Use Case:** Recommended for $ t < 10,000 $ or when checking for high accuracy.
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
         * * **Precision:** High (depends on Taylor series depth).
         * * **Use Case:** Only efficient when computing **blocks** of zeros for extremely large $ t $ (e.g., $ t > 10^{10} $).
         */
        OdlyzkoSchonhage
    };

    /**
     * @class HardyZ
     * @brief A static class for computing the Hardy Z-function on the critical line.
     * * The Hardy Z-function is defined as:
     * $$
     * Z(t) = e^{i \theta(t)} \zeta\left(\frac{1}{2} + it\right)
     * $$
     * where $ \theta(t) $ is the Riemann-Siegel theta function and 
     * $ \zeta(s) $ is the Riemann zeta function.
     * * By definition, $ Z(t) $ is real-valued for real $ t $.
     * $ |Z(t)| = |\zeta(1/2 + it)| $, making it useful for studying the zeros of Zeta.
     */
    class HardyZ {
    public:
        // =============================================================
        // Single Point Computation (EM and RS)
        // =============================================================
        /**
         * @brief Computes the value of Z(t).
         * @param t The imaginary component of the argument $ s = \frac{1}{2} + it $.
         * @param method The algorithm to use (default: EulerMaclaurin).
         * @return The real value $ Z(t) $.
         */
        static double compute(double t, Method method = Method::EulerMaclaurin);

        // =============================================================
        // Block Computation (Odlyzko-Schönhage)
        // =============================================================
        /**
         * @brief Computes a range of Z values efficiently.
         * This is the primary entry point for the Odlyzko-Schönhage method.
         * It computes $ Z(t) $ for $ t \in [start, start + length] $.
         * @param start_t The starting height on the critical line.
         * @param length The length of the interval to compute.
         * @param points The number of sampling points within the interval.
         * @param method Algorithm (default: OdlyzkoSchonhage for efficiency).
         * @return A vector containing the Z values.
         */
        static std::vector<double> computeBlock(double start_t, double length, int points, Method method = Method::OdlyzkoSchonhage);   


    private:
        /**
         * @brief Helper for computeEM to calculate the complex Zeta value.
         * Computes $ \zeta(s) $ via:
         * $$
         * \zeta(s) \approx \sum_{n=1}^{N-1} n^{-s} + \frac{N^{1-s}}{s-1} + \frac{1}{2}N^{-s} + \sum_{k=1}^{m} \frac{B_{2k}}{(2k)!} f^{(2k-1)}(N)
         * $$
         * @param s The complex argument $ \frac{1}{2} + it $.
         * @param N The summation cutoff limit.
         */
        static std::complex<double> zetaEM(std::complex<double> s, int N);

        /**
         * @brief Computes Z(t) using the Euler-Maclaurin summation for $ \zeta(s) $.
         * Computes $ \zeta(s) $ via zetaEM 
         * and rotates the result by $ e^{i\theta(t)} $.
         */
        static double computeEM(double t);

        /**
         * @brief Computes Z(t) using the Riemann-Siegel Main Sum.
         * Note that it does not use $ \zeta(s) $ as the Euler-Maclaurin method did.
         * Uses the approximation:
         * $$
         * Z(t) \approx 2 \sum_{n=1}^{\lfloor \sqrt{t/2\pi} \rfloor} \frac{\cos(\theta(t) - t \ln n)}{\sqrt{n}}
         * $$
         */
        static double computeRS(double t);
    
        /**
         * @brief Implementation of the rational function expansion / FFT method.
         * * Uses the Taylor expansion of $ n^{-i\delta} $ to evaluate sums quickly.
         */
        static std::vector<double> computeOS(double start_t, double length, int points);
    };
}

#endif 