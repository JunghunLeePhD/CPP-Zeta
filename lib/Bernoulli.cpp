#include "Bernoulli.h"
#include <cmath>
#include <numeric> 
#include <vector>

namespace Zeta {

    std::vector<double> Bernoulli::cache = { 1.0 };

    double Bernoulli::nCr(int n, int k) {
        if (k < 0 || k > n) return 0.0;
        if (k == 0 || k == n) return 1.0;
        if (k > n / 2) k = n - k;
        
        std::vector<int> steps(k);
        std::iota(steps.begin(), steps.end(), 1); 

        return std::accumulate(
            steps.begin(), steps.end(),
            1.0, 
            [n](double current_res, int i) {
                return current_res * (n - i + 1) / i;
            }
        );
    }

    double Bernoulli::get(int n) {
        if (n < cache.size()) {
            return cache[n];
        }

        for (int m = cache.size(); m <= n; ++m) {
            
            std::vector<int> k_range(m);
            std::iota(k_range.begin(), k_range.end(), 0);

            double sum = std::accumulate(
                k_range.begin(), k_range.end(),
                0.0,
                [m](double current_sum, int k) { 
                    double term = nCr(m + 1, k) * cache[k];
                    return current_sum + term; 
                }
            );
            double val = -1.0 / (m + 1.0) * sum;
            cache.push_back(val);
        }

        return cache[n];
    }
}