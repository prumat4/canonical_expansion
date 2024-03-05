#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <bitset>
#include <map>
#include <algorithm>

std::map<uint64_t, std::vector<uint64_t>> result_map;

void precompute_result_map() {
    std::vector<uint64_t> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101};

    for (uint64_t d : primes) {
        result_map[d] = {1};
        while (std::count(result_map.at(d).begin(), result_map.at(d).end(), result_map.at(d).back()) < 2) {
            result_map.at(d).push_back((result_map.at(d).back() * 2) % d);
        }

        result_map[d].pop_back();
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

uint64_t method_of_trial_divisions(uint64_t number) {
    std::bitset<64> binary_number(number);

    if (!binary_number.test(0)) {
        return 2;
    }

    for (auto it = result_map.begin(); it != result_map.end(); ++it) {
        uint64_t d = it->first;
        uint64_t sum = 0;

        // problem lies here
        for (size_t i = 63; i > 0; i--)
            sum += binary_number.test(i) * result_map.at(d).at(i % result_map.at(d).size());
            sum %= d;

        if (sum == 0)
            return d;
    }

    return 1;
}

int main() {
    precompute_result_map();
    // print_result_map();
    std::cout << method_of_trial_divisions(67443) << std::endl;

    return 0;
}