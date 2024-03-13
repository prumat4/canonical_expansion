#pragma once

#include <cstdint>
#include <cstdlib>
#include <vector>

int64_t evaluateJacobiSymbol(uint64_t base, uint64_t value);
bool solovey_strassen_test(const long number, const std::size_t iterationsNumber = 100);
std::vector<uint64_t> generateFactorBase(const uint64_t number, const int limit);
std::vector<uint64_t> continuedFraction(const uint64_t number);
std::vector<uint64_t> findSmoothValues(const uint64_t number);