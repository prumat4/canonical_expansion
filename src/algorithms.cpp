#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <gmpxx.h>

int64_t evaluateJacobiSymbol(uint64_t base, uint64_t value)
{
    if(std::gcd(value, base) != 1)
        return 0;

    if(value <= 0 || value % 2 != 1)
        return 0;
    
    base %= value;
    int64_t t = 1;
    int64_t r;

    while (base != 0) 
    {
        while (base % 2 == 0) 
        {
            base /= 2;
            r = value % 8;
            if (r == 3 || r == 5) 
                t = -t;
        }

        r = value;
        value = base;
        base = r;

        if (base % 4 == 3 && value % 4 == 3)
            t = -t;

        base %= value;
    }

    if(value == 1)
        return t;
    else
        return 0;
}

bool solovey_strassen_test(const uint64_t number, const std::size_t iterationsNumber = 100)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> distribution(1, number - 1);

    for(std::size_t i = 0; i < iterationsNumber; ++i)
    {
        auto randomValue = distribution(gen);
        uint64_t possibleDivider = std::gcd(number, randomValue);

        if(possibleDivider != 1)
            return false;

        int64_t jacobiSymbol = evaluateJacobiSymbol(randomValue, number);
        uint64_t power = std::pow(randomValue, (number - 1) / 2);

        while(jacobiSymbol < 0)
            jacobiSymbol += number;
        power %= number;

        if(jacobiSymbol != power)
        {
            return false;
        }
    }

    return true;
}

std::vector<uint64_t> generateFactorBase(const uint64_t number, const int limit)
{
    std::size_t k = 100; // hardcode
    std::vector<bool> isPrime(limit + 1, true);
    std::vector<uint64_t> factorBase;

    std::random_device rd;
    std::mt19937 gen(rd());

    const long double limitor = std::exp(std::pow(std::log(number) *
                                std::log(std::log(number)), 1.0/2.0));
    const double alpha = 1.0 / std::sqrt(2);
    std::uniform_int_distribution<uint64_t> distribution(1, std::pow(limitor, alpha));
    std::cout << "***" << std::pow(limitor, alpha) << std::endl;

    for(std::size_t i = 0; i < limit && factorBase.size() < k; ++i)
    {
        uint64_t possiblePrimeNumber = distribution(gen);
        std::cout << "possiblePrimer: " << possiblePrimeNumber << std::endl;

        if(solovey_strassen_test(possiblePrimeNumber))
        {
            if(std::find(factorBase.begin(), factorBase.end(), possiblePrimeNumber) == factorBase.end())
                factorBase.push_back(possiblePrimeNumber);
        }
    }

    return factorBase;
}