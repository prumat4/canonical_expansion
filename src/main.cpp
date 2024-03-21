#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <bitset>
#include <map>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <random>
#include "algorithms.h"

std::map<uint64_t, std::vector<uint64_t>> result_map;

void precompute_result_map() {
    std::vector<uint64_t> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

    for (uint64_t d : primes) {
        result_map[d] = {1};
        while (std::count(result_map.at(d).begin(), result_map.at(d).end(), result_map.at(d).back()) < 2) {
            result_map.at(d).push_back((result_map.at(d).back() * 2) % d);
        }
        result_map.at(d).pop_back();
    }
}

void print_result_map() {
    for (const auto& pair : result_map) {
        std::cout << pair.first << ": ";
        for (uint64_t num : pair.second) {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }
}

uint32_t method_of_trial_divisions(const uint64_t number) {
    if (number < 2) return number;

    if ((number & 1) == 0) return 2;

    for (const auto& it : result_map) {
        uint64_t d = it.first;
        uint64_t sum = 0;

        for (size_t i = 0; i < 64; ++i) {
            bool bitSet = (number >> i) & 1;
            sum = (sum + (bitSet ? result_map[d][i % result_map[d].size()] : 0)) % d;
        }

        if (sum == 0) return d;
    }

    return 1;
}

int rho_function(int x, int modulus) {
    return (x * x + 1) % modulus;
}

uint64_t rho_pollard(uint64_t number) {
    std::vector<uint64_t> x_vector;
    std::vector<uint64_t> y_vector;

    x_vector.push_back(rho_function(2, number));
    y_vector.push_back(rho_function(x_vector.at(0), number));

    uint64_t i = 0;
    while (std::count(x_vector.begin(), x_vector.end(), x_vector.at(i)) != 2) {
        x_vector.push_back(rho_function(x_vector.at(i), number));
        y_vector.push_back(rho_function(rho_function(y_vector.at(i), number), number));
        i++;

        if (x_vector.at(i) - y_vector.at(i) == 0) 
            return 1;
    
        uint64_t divisor = std::gcd(number, std::abs(static_cast<int64_t>(y_vector.at(i) - x_vector.at(i))) % number);
        if (divisor != 1) 
            return divisor;
    }

    return 1;
}

int main(int argc, char** argv) {
    if (argc > 1) {
        int64_t number = std::stoll(argv[1]);
        std::cout << number << std::endl;
    } else {
        std::cerr << "Error: No number provided.\n"
                  << "Usage: " << argv[0] << " <number>\n\n"
                  << "Please provide a number as an argument. This number will be used for [briefly describe what the number is for].\n"
                  << "Example:\n"
                  << "  " << argv[0] << " 12345\n\n"
                  << "In this example, '12345' is the number that will be used to [briefly describe what the program does with the number].\n";
        // should think about this later :)
        // return 1;
    }
    
    precompute_result_map();
    std::vector<uint64_t> test_values = { 31, 1021514194991569, 499664789704823, 3009182572376191, 666666, 99999999, 123456789, 872, 816258 };

    for(const auto value : test_values) {
        uint32_t divisor = method_of_trial_divisions(value);
        if (divisor == 1) {
            std::cout << value << " is prime or not divisible by the primes in the map.\n";
        } else {
            std::cout << std::fixed << std::setprecision(4);
            std::cout << value << " is divisible by " << divisor << ", " << value << " / " << divisor << " = " << float(value) / float(divisor) << std::endl;
        }
        divisor = rho_pollard(value);
        if (divisor == 1) {
            std::cout << value << " is prime" << std::endl;
        } else {
            std::cout << value << " is divisible by " << divisor << ", " << value << " / " << divisor << " = " << value / divisor << std::endl;
        }
    }

    // std::cout << evaluateJacobiSymbol(2, 3) << std::endl;
    // auto test = generateFactorBase(1237812319023, 1000000);
    // for(const auto& it : test)
    // {
    //     std::cout << it << " " << std::endl;
    // }
    // std::cout << "solovey: " << solovey_strassen_test(31) << std::endl;
    
    // auto k = methodBrilhartMorrison(72365723);
    // std::cout << k;
    
    return 0;
}