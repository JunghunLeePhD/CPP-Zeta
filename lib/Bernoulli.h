#pragma once

#include <vector>
#include <concepts>

namespace Zeta {

    /**
     * @brief Retrieves the n-th Bernoulli number B_n.
     * Uses a local static cache to store previously computed values (Memoization).
     * Formula:
     * $$ B_m = \frac{-1}{m+1} \sum_{k=0}^{m-1} \binom{m+1}{k} B_k $$
     * @tparam T Floating point type (float, double, long double).
     * @param n The index (must be >= 0).
     * @return The Bernoulli number.
     */
    template <std::floating_point T>
    [[nodiscard]]
    T bernoulli(int n);

} 

#include "Bernoulli.tpp"