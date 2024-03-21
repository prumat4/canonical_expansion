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

int main(int argc, char** argv) {
    if (argc > 1) {
        try {
            uint64_t number = std::stoull(argv[1]);
            std::cout << "Processing number: " << number << std::endl;

            precomputeResultMap(); // Ensure any necessary setup is done before computations

            // Call the function to find the canonical expansion and print the factors
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