#include <cmath>
#include <numbers> 

namespace Zeta {

    template <std::floating_point T>
    constexpr T theta(T t) noexcept {
        constexpr T PI = std::numbers::pi_v<T>;

        if (std::abs(t) < T{1e-9}) return T{0}; 

        const T half_t = t / T{2};
        
        const T term_log    = half_t * std::log(t / (T{2} * PI));
        const T term_linear = -half_t;
        const T term_const  = -PI / T{8};
        const T term_corr1  = T{1} / (T{48} * t);
        
        const T t3          = t * t * t;
        const T term_corr2  = T{7} / (T{5760} * t3);

        return term_log + term_linear + term_const + term_corr1 + term_corr2;
    }

} 