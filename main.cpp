#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <bitset>
#include <map>
#include <algorithm>
#include <numeric>

std::map<uint64_t, std::vector<uint64_t>> result_map;

void precompute_result_map() {
    std::vector<uint64_t> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101};

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
    std::bitset<64> binary_number(number);

    if (binary_number.test(binary_number.size() - 1))
        return 2;

       for (auto it = result_map.begin(); it != result_map.end(); ++it) {
        uint64_t d = it->first;
        uint64_t sum = 0;

        for (size_t i = 0; i < binary_number.size(); i++) {
            // m?
            sum += (binary_number[binary_number.size() - i - 1]) * result_map.at(d).at(i % result_map.at(d).size());
        }

        if (sum % d == 0)
            return d;
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
        

        uint64_t divisor = std::gcd(number, abs(y_vector.at(i) - x_vector.at(i)) % number);
        if (divisor != 1) 
            return divisor;
    }

    return 1;
}

int main() {
    precompute_result_map();
    std::vector<uint64_t> test_values = { 1021514194991569, 499664789704823, 3009182572376191, 666666, 99999999, 123456789 };

    for(const auto value : test_values) {
        std::cout << method_of_trial_divisions(value) << std::endl;

        uint64_t divisor = rho_pollard(value);
        if (divisor == 1) {
            std::cout << value << " is prime" << std::endl;
        } else {
            std::cout << value << " is divisible by " << divisor << ", " << value << " / " << divisor << " = " << value / divisor << std::endl;
        }
    }

    return 0;
}