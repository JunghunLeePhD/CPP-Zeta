#ifndef BERNOULLI_H
#define BERNOULLI_H

#include <vector>

namespace Zeta {

    class Bernoulli {
    public:
        /**
         * @brief Retrieves the n-th Bernoulli number B_n.
         * * Uses a cache to store previously computed values.
         * * Computed using the recursive formula: 
         * $ \sum_{k=0}^n \binom{n+1}{k} B_k = 0 $
         * That is,
         * $ B_m = -1/(m+1) * Sum(k=0 to m-1) [ binom(m+1, k) * B_k ] $
         * @param n The index (must be >= 0).
         * @return The Bernoulli number as a double.
         */
        static double get(int n);

    private:
        static std::vector<double> cache;
        
        static double nCr(int n, int k);
    };

}

#endif