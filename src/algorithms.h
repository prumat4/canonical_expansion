#pragma once

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <map>

std::vector<std::vector<int64_t>> readMatrix(int rows, int cols);
void printSolution(const std::vector<int64_t>& solution);

int64_t evaluateJacobiSymbol(uint64_t base, uint64_t value);
bool solovey_strassen_test(const long number, const std::size_t iterationsNumber = 100);
std::vector<int64_t> generateFactorBase(const uint64_t number, const double alpha);
std::vector<int64_t> continuedFraction(const uint64_t number, const std::size_t k);
std::vector<int64_t> findSmoothValues(const uint64_t number, std::vector<int64_t> factorBase, const std::size_t k);
std::pair<std::vector<int64_t>, bool> factorizeSmoothValue(const int64_t smoothValue, const int64_t number, const std::vector<int64_t>& factorBase);

int64_t methodBrilhartMorrison(const uint64_t number);