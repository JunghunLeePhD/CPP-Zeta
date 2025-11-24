#include <iostream>
#include <algorithm>
#include <vector>
#include "Theta.h"
#include "HardyZ.h"

int main() {
    int size = 500;
    double start_value = 0.0;
    double step = 0.1f;

    std::vector<double> v(size);
    std::generate(v.begin(), v.end(), [current = start_value, &step]() mutable{
        double value = current;
        current += step;

        double hardyZ = Zeta::HardyZ::compute(value, Zeta::Method::EulerMaclaurin);
        std::cout << hardyZ << std::endl;

        return value;
    });

    return 0;
}