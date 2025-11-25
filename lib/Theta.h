#pragma once

#include <concepts>

namespace Zeta {

    /**
     * @brief Computes the Riemann-Siegel Theta function $ \theta(t) $.
     * Defined as:
     * $$ \theta(t) = \text{Im} \left[ \ln \Gamma\left(\frac{1}{4} + i\frac{t}{2}\right) \right] - \frac{t}{2} \ln \pi $$
     * Uses Stirling's asymptotic expansion for large $ t $.
     * @tparam T Floating point type (float, double, long double).
     * @param t The imaginary part of the argument.
     * @return The phase angle value.
     */
    template <std::floating_point T>
    [[nodiscard]] 
    constexpr T theta(T t) noexcept;

} 

#include "Theta.tpp"