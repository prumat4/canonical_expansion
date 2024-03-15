#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
// #include <gmp.h>
#include <gmpxx.h>
#include <stdio.h>

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

bool solovey_strassen_test(const long number, const std::size_t iterationsNumber = 100) 
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<long> distribution(3, number - 1);

    mpz_class mpzNumber(number);
    for (std::size_t i = 0; i < iterationsNumber; ++i) {
        mpz_class randomValue = distribution(gen);
        mpz_class possibleDivider;
        mpz_gcd(possibleDivider.get_mpz_t(), mpzNumber.get_mpz_t(), randomValue.get_mpz_t());

        // Check if gcd(randomValue, number) = 1
        if (mpz_cmp_ui(possibleDivider.get_mpz_t(), 1) != 0)
            return false;

        if (mpz_odd_p(possibleDivider.get_mpz_t()) == 0)
            return false;

        int64_t jacobiSymbol = mpz_jacobi(randomValue.get_mpz_t(), mpzNumber.get_mpz_t());
        mpz_class power;
        mpz_powm_ui(power.get_mpz_t(), randomValue.get_mpz_t(), (number - 1) / 2, mpzNumber.get_mpz_t());

        while (jacobiSymbol < 0)
            jacobiSymbol += number;
        power %= mpzNumber;

        if (jacobiSymbol != power.get_si())
            return false;
    }

    return true;
}

std::vector<int64_t> sieveOfEratosthenes(const double limit) 
{
    std::vector<bool> isPrime(limit + 1, true);
    isPrime[0] = isPrime[1] = false;

    for (int i = 2; i * i <= limit; ++i) 
    {
        if (isPrime[i]) 
        {
            for (int j = i * i; j <= limit; j += i) 
            {
                isPrime[j] = false;
            }
        }
    }

    std::vector<int64_t> factorBase;
    factorBase.push_back(-1);

    for (int i = 2; i <= limit; ++i) 
    {
        if (isPrime[i]) 
        {
            mpz_class mpzValue = static_cast<int>(limit);
            mpz_class mpzI = i;
            if(mpz_legendre(mpzValue.get_mpz_t(), mpzI.get_mpz_t()) == 1)
                factorBase.push_back(i);
        }
    }

    return factorBase;
}

std::vector<uint64_t> generateFactorBase(const uint64_t number, const int limit)
{
    
    // std::size_t k = 100; // hardcode
    // std::vector<bool> isPrime(limit + 1, true);
    // std::vector<uint64_t> factorBase;

    // std::random_device rd;
    // std::mt19937 gen(rd());
    int64_t newNumber = 17873;
    const long double limitor = std::exp(std::pow(std::log(newNumber) *
                                std::log(std::log(newNumber)), 1.0/2.0));
    const double alpha = 1.0 / std::sqrt(2);

    std::vector<int64_t> factorBase = sieveOfEratosthenes(std::pow(limitor, alpha));
    for(const auto& it : factorBase)
    {
        std::cout << it << " ";
    }
    // std::uniform_int_distribution<uint64_t> distribution(1, std::pow(limitor, alpha));

    // for(std::size_t i = 0; i < limit && factorBase.size() < k; ++i)
    // {
    //     int64_t possiblePrimeNumber = distribution(gen);

    //     if(solovey_strassen_test(possiblePrimeNumber))
    //     {
    //         if(std::find(factorBase.begin(), factorBase.end(), possiblePrimeNumber) == factorBase.end())
    //         {
    //             mpz_class mpzNumber = static_cast<int>(number);
    //             mpz_class mpzPossiblePrimeNumber = static_cast<int>(possiblePrimeNumber);
    //             if(mpz_legendre(mpzNumber.get_mpz_t(), mpzPossiblePrimeNumber.get_mpz_t()) == 1)
    //                 factorBase.push_back(possiblePrimeNumber);
    //         }
    //     }
    // }

    // return factorBase;

    return std::vector<uint64_t>();
}

std::vector<uint64_t> continuedFraction(const uint64_t number)
{
    std::vector<uint64_t> values;

    std::size_t k = 100; // hardcode
    double squareNumber = std::sqrt(number);

    uint64_t v_i = 1;
    uint64_t alpha_i = squareNumber;
    uint64_t a_i = std::floor(alpha_i);
    uint64_t u_i = a_i;
    
    for(std::size_t i = 0; i < k + 1; ++i)
    {
        values.push_back(a_i);

        v_i = (number - (u_i * u_i)) / v_i;
        alpha_i = (squareNumber + u_i) / v_i;
        a_i = std::floor(alpha_i);
        u_i = a_i * v_i - u_i;
    }
    values.push_back(a_i);

    return values;
}

std::vector<uint64_t> findSmoothValues(const uint64_t number)
{
    std::size_t k = 100; // hardcode
    std::vector<uint64_t> smoothValues;

    auto continuedFractionValues = continuedFraction(number);

    int64_t b_i = 1;
    int64_t b_ii = 0;

    for(std::size_t i = 0; i < k + 1; ++i)
    {
        uint64_t result = continuedFractionValues[i] * b_i + b_ii;
        result = (result * result) % number;

        smoothValues.push_back(result);

        b_ii = b_i;
        b_i = result;
    }

    for(const auto& it : smoothValues)
    {
        std::cout << it << " ";
    }
    std::cout << std::endl << "size: " << smoothValues.size();

    return smoothValues;
}