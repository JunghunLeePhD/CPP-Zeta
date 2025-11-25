#include <print>
#include "Theta.h"

int main() {
    std::println("--- Riemann-Siegel Theta Function Test (C++26) ---");

    // 1. Standard Runtime Test (Double)
    double t_double = 14.13472514173469; 
    double res_double = Zeta::theta(t_double);
    
    std::println("Type: Double");
    std::println("t = {}", t_double);
    std::println("Result = {}\n", res_double);

    // 2. Float Test (Template instantiation)
    float t_float = 21.022040f; 
    float res_float = Zeta::theta(t_float);

    std::println("Type: Float");
    std::println("t = {}", t_float);
    std::println("Result = {}\n", res_float);

    // 3. Constexpr Test (Compile-Time Calculation)
    const double t_const = 10.0; 
    const double val_const = Zeta::theta(t_const);

    std::println("Type: Constexpr (Computed at Compile-Time)");
    std::println("t = {}", t_const);
    std::println("Result = {}", val_const); 

    return 0;
    return 0;
}