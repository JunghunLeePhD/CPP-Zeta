#include <iostream>
#include <algorithm>
#include <vector>
#include "Theta.h"
#include "HardyZ.h"

int main() {
    int size = 500;
    double start_value = 100000.0;
    double step = 0.1f;

    std::vector<double> v(size);
    std::generate(v.begin(), v.end(), [current = start_value, &step]() mutable{
        double value = current;
        current += step;

        double hardyZ_EM = Zeta::HardyZ::compute(value, Zeta::Method::EulerMaclaurin);
        double hardyZ_RS = Zeta::HardyZ::compute(value, Zeta::Method::RiemannSiegel);

        std::cout << hardyZ_EM << " " << hardyZ_RS << std::endl;

        return value;
    });

    return 0;
}