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
#include <chrono>

void testRhoPollardPerformance() {
    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<uint64_t> dist(std::llround(std::pow(2,10)), std::llround(std::pow(2,20)));

    std::vector<double> executionTimes;
    for (int i = 0; i < 100; ++i) {
        uint64_t number = dist(e2);

        auto start = std::chrono::high_resolution_clock::now();
        auto result = rhoPollard(number);
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> elapsed = end - start;
        executionTimes.push_back(elapsed.count());
        std::cout << "Result for " << number << ": " << result << ", Time: " << elapsed.count() << " ms\n";
    }

    double averageTime = std::accumulate(executionTimes.begin(), executionTimes.end(), 0.0) / executionTimes.size();
    std::cout << "Average execution time: " << averageTime << " ms\n";
}


int main(int argc, char** argv) {
    if (argc > 1) {
        try {
            uint64_t number = std::stoull(argv[1]);
            std::cout << "Processing number: " << number << std::endl;
            precomputeResultMap();

            testRhoPollardPerformance();

            std::vector<uint64_t> factors = findCanonicalExpansion(number);
            std::cout << "Canonical expansion: ";
            for (const auto& factor : factors) {
                std::cout << factor << " ";
            }
            std::cout << std::endl;

        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid number format.\n";
            return 1;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Number out of range.\n";
            return 1;
        }
    } else {
        std::cerr << "Error: No number provided.\n"
                  << "Usage: " << argv[0] << " <number>\n\n"
                  << "Please provide a number as an argument. This number will be used to find its canonical expansion.\n"
                  << "Example:\n"
                  << "  " << argv[0] << " 12345\n\n"
                  << "In this example, '12345' is the number that will be analyzed to find and print its canonical expansion.\n";
        return 1;
    }
    
    return 0;
}