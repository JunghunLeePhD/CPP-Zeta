#include <print>
#include <algorithm> 
#include <ranges>
#include "HardyZ.h"

int main(){
    std::println("==========================================");
    std::println("   C++26 Zeta Function Library Test");
    std::println("==========================================\n");

    // Known 1st Zero: t = 14.13472514173469...
    double t_zero = 14.13472514173469;

    std::println("    Target t = {:.14f}", t_zero);

    // Method A: Euler-Maclaurin (High Precision)
    double z_em = Zeta::Hardy::compute<double>(t_zero, Zeta::Method::EulerMaclaurin);
    std::println("    [EM] Z(t) = {: .15f} (Expected ~0.0)", z_em);

    // Method B: Riemann-Siegel (Approximation)
    double z_rs = Zeta::Hardy::compute<double>(t_zero, Zeta::Method::RiemannSiegel);
    std::println("    [RS] Z(t) = {: .15f} (Expected ~0.0)", z_rs);
    std::println("");

    // Method C: Odlyzko-Schönhage 
    double start_t = t_zero - 3;
    double length = 7.0;  // Compute interval [100.0, 105.0]
    int points = 8;       // 8 sample points

    std::println("    Range: [{}, {}], Points: {}", start_t, start_t + length, points);

    auto results = Zeta::Hardy::computeBlock<double>(
        start_t, length, points, Zeta::Method::OdlyzkoSchonhage
    );

    int idx = 0;
    double step = length / (points - 1);
    
    for (auto [idx, val] : results | std::views::enumerate) {
        double current_t = start_t + (static_cast<double>(idx) * step);
        std::println("    [OS] Z(t) = {: .15f}  <-  t = {:.14f}", val, current_t);
    }

    return 0;
} 