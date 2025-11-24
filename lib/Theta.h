#ifndef THETA_H
#define THETA_H

namespace Zeta {

    /**
     * @class Theta
     * @brief Helper class to compute the Riemann-Siegel Theta function.
     */
    class Theta {
    public:
        /**
         * @brief Computes the Riemann-Siegel Theta function $ \theta(t) $.
         * * Defined as:
         * $$
         * \theta(t) = \text{Im} \left[ \ln \Gamma\left(\frac{1}{4} + i\frac{t}{2}\right) \right] - \frac{t}{2} \ln \pi.
         * $$
         * This implementation uses Stirling's asymptotic expansion for large $ t $:
         * $$
         * \theta(t) \approx \frac{t}{2} \ln\left(\frac{t}{2\pi}\right) - \frac{t}{2} - \frac{\pi}{8} + \frac{1}{48t} + \dots.
         * $$
         * @param t The imaginary part of the argument.
         * @return The phase angle value.
         */
        static double value(double t);
    };

}

#endif