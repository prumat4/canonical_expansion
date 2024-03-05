#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <bitset>
#include <map>
#include <algorithm>

std::map<int, std::vector<int>> result_map;

void precompute_result_map() {
    std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101};

    for (int d : primes) {
        result_map[d] = {1};
        while (std::count(result_map.at(d).begin(), result_map.at(d).end(), result_map.at(d).back()) < 2) {
            result_map.at(d).push_back((result_map.at(d).back() * 2) % d);
        }

        result_map[d].pop_back();
    }
}

int main() {
    precompute_result_map();

    for (const auto& pair : result_map) {
        std::cout << pair.first << ": ";
        for (int num : pair.second) {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}