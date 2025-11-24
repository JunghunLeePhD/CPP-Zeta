#include <iostream>
#include <algorithm>
#include <vector>
#include "Theta.h"

int main() {
    int size = 10;
    double start_value = 0.0;
    double step = 0.1f;

    std::vector<double> v(size);
    std::generate(v.begin(), v.end(), [current = start_value, &step]() mutable{
        double value = current;
        current += step;

        double theta = Zeta::Theta::value(value);
        std::cout << theta << std::endl;

        return value;
    });

    return 0;
}