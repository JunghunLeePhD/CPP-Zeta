#include <print>
#include <algorithm> 
#include <ranges>
#include "Bernoulli.h"

int main(){
    std::println("--- Bernoulli Numbers (Memoized) ---");

    std::ranges::for_each(
        std::views::iota(0, 13), 
        [](int i) {
            double b = Zeta::bernoulli<double>(i);
            std::println("B_{:<2} = {: .10f}", i, b);
        }
    );

    std::println("--- Template Instantiation Check ---");
    
    float b_float = Zeta::bernoulli<float>(2); 
    std::println("Float B_2  = {}", b_float);
} 