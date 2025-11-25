#include <vector>
#include <cmath>
#include <numeric>    
#include <ranges>    
#include <algorithm> 

namespace Zeta {

    /**
     * @brief Computes the Binomial Coefficient $\binom{n}{k}$ (n choose k).
     * **Definition:**
     * $$ \binom{n}{k} = \frac{n!}{k!(n-k)!} $$
     * **Implementation:**
     * Uses the multiplicative formula to avoid calculating large factorials directly:
     * $$ \binom{n}{k} = \prod_{i=1}^k \frac{n - i + 1}{i} $$
     * @note Marked `constexpr`. If $n$ and $k$ are compile-time constants, 
     * this value is computed by the compiler (Zero Runtime Cost).
     * @param n Total number of items.
     * @param k Number of items to choose.
     * @return The binomial coefficient as a double.
     */
    constexpr double nCr_impl(int n, int k) {
        if (k < 0 || k > n) return 0.0;
        if (k == 0 || k == n) return 1.0;

        if (k > n / 2) k = n - k;

        auto range = std::views::iota(1, k + 1);

        return std::ranges::fold_left(
            range, 
            1.0, 
            [n](double acc, int i) {
                return acc * (n - i + 1) / i;
            }
        );
    }

    template <std::floating_point T>
    T bernoulli(int n) {
        static std::vector<T> cache = { T{1.0} }; // B_0 = 1

        if (n < static_cast<int>(cache.size())) return cache[n];

        for (int m = cache.size(); m <= n; ++m) {
            auto k_range = std::views::iota(0, m);

            T sum = std::ranges::fold_left(
                k_range,
                T{0},
                [m](T acc, int k) {
                    T combinations = static_cast<T>(nCr_impl(m + 1, k));
                    return acc + (combinations * cache[k]);
                }
            );

            T b_m = (T{-1} / static_cast<T>(m + 1)) * sum; // B_m = -1/(m+1) * sum
            cache.push_back(b_m);
        }

        return cache[n];
    }

} 