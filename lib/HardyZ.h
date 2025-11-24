#ifndef HARDYZ_H
#define HARDYZ_H

#include <complex>

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
        /**
         * @brief Computes the value of Z(t).
         * @param t The imaginary component of the argument $ s = \frac{1}{2} + it $.
         * @param method The algorithm to use (default: EulerMaclaurin).
         * @return The real value $ Z(t) $.
         */
        static double compute(double t, Method method = Method::EulerMaclaurin);

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
    };
}

#endif 