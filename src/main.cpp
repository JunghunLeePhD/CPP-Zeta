#include <iostream>
#include <iomanip>
#include "HardyZ.h"

int main() {
    double t_start = 1000.0;
    double length = 5.0;
    int points = 6; // Evaluate 1000, 1001, ... 1005

    std::cout << "Comparing RS (Single) vs Odlyzko (Block)" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    auto block_results = Zeta::HardyZ::computeBlock(
        t_start, length, points, Zeta::Method::OdlyzkoSchonhage
    );

    double step = length / (points - 1);

    for(int i = 0; i < points; ++i) {
        double t = t_start + i * step;
        double val_os = block_results[i];
        double val_rs = Zeta::HardyZ::compute(t, Zeta::Method::RiemannSiegel);

        std::cout << "t=" << std::fixed << std::setprecision(1) << t 
                  << " | OS: " << std::setprecision(6) << val_os 
                  << " | RS: " << val_rs 
                  << " | Diff: " << std::scientific << std::abs(val_os - val_rs) 
                  << std::endl;
    }

    return 0;
}