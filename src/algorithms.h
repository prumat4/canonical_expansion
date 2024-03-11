#pragma once

#include <cstdint>
#include <cstdlib>
#include <vector>

int64_t evaluateJacobiSymbol(uint64_t base, uint64_t value);
bool solovey_strassen_test(const uint64_t number, const std::size_t iterationsNumber = 100);
std::vector<uint64_t> generateFactorBase(const uint64_t number, const int limit);