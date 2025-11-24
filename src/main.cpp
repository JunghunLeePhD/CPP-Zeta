#include <iostream>
#include <algorithm>
#include <vector>
#include "Theta.h"
#include "HardyZ.h"
#include "Bernoulli.h"


int main() {
    int size = 10;
    int start_value = 0;
    int step = 1;

    std::vector<int> v(size);
    std::generate(v.begin(), v.end(), [current = start_value, &step]() mutable{
        int value = current;
        current += step;
        double bernoulli = Zeta::Bernoulli::get(value);

        std::cout << value << " " << bernoulli << std::endl;

        return value;
    });

    return 0;
}