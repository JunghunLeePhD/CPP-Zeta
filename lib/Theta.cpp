#include "Theta.h"
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Zeta {

    double Theta::value(double t) {
        if (std::abs(t) < 1e-9) return 0.0; 

        double half_t = t / 2.0;
        double term_log = half_t * std::log(t / (2.0 * M_PI));
        double term_linear = -half_t;
        double term_const = -M_PI / 8.0;
        double term_corr1 = 1.0 / (48.0 * t);
        double t3 = t * t * t;
        double term_corr2 = 7.0 / (5760.0 * t3);

        return term_log + term_linear + term_const + term_corr1 + term_corr2;
    }

}