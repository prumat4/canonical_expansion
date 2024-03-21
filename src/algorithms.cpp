#include <map>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <gmpxx.h>
#include <stdio.h>
#include <chrono>
#include <utility>

std::map<uint64_t, std::vector<uint64_t>> result_map;

void precomputeResultMap() {
    std::vector<uint64_t> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    for (uint64_t d : primes) {
        result_map[d] = {1};
       
        while (std::count(result_map.at(d).begin(), result_map.at(d).end(), result_map.at(d).back()) < 2) {
            result_map.at(d).push_back((result_map.at(d).back() * 2) % d);
        }
        result_map.at(d).pop_back();
    }
}

void printResultMap() {
    for (const auto& pair : result_map) {
        std::cout << pair.first << ": ";
        for (uint64_t num : pair.second) {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }
}

uint32_t methodOfTrialDivisions(const uint64_t number) {
    if (number < 2) return number;
   
    if ((number & 1) == 0) return 2;
   
    for (const auto& it : result_map) {
        uint64_t d = it.first;
        uint64_t sum = 0;
        for (size_t i = 0; i < 64; ++i) {
            bool bitSet = (number >> i) & 1;
            sum = (sum + (bitSet ? result_map[d][i % result_map[d].size()] : 0)) % d;
        }
        if (sum == 0) return d;
    }
    return 1;
}

int rhoFunction(int x, int modulus) {
    return (x * x + 1) % modulus;
}

uint64_t rhoPollard(uint64_t number) {
    std::vector<uint64_t> x_vector;
    std::vector<uint64_t> y_vector;
  
    x_vector.push_back(rhoFunction(2, number));
    y_vector.push_back(rhoFunction(x_vector.at(0), number));
  
    uint64_t i = 0;
  
    while (std::count(x_vector.begin(), x_vector.end(), x_vector.at(i)) != 2) {
        x_vector.push_back(rhoFunction(x_vector.at(i), number));
        y_vector.push_back(rhoFunction(rhoFunction(y_vector.at(i), number), number));
        i++;
      
        if (x_vector.at(i) - y_vector.at(i) == 0) 
            return 1;
  
        uint64_t divisor = std::gcd(number, std::abs(static_cast<int64_t>(y_vector.at(i) - x_vector.at(i))) % number);
  
        if (divisor != 1) 
            return divisor;
    }
    return 1;
}

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

std::vector<int64_t> sieveOfEratosthenes(const uint64_t number, const double limit) 
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

    mpz_class mpzValue = static_cast<int>(number);

    std::cout << "limit: " << limit << std::endl;
    for (int i = 2; i <= limit; ++i) 
    {
        if (isPrime[i]) 
        {
            mpz_class mpzI = i;
            if(mpz_legendre(mpzValue.get_mpz_t(), mpzI.get_mpz_t()) == 1)
                factorBase.push_back(i);
        }
    }

    return factorBase;
}

std::vector<int64_t> generateFactorBase(const uint64_t number, const double alpha)
{
    const long double limitor = std::exp(std::pow(std::log(number) *
                                std::log(std::log(number)), 1.0/2.0));

    return sieveOfEratosthenes(number, std::pow(limitor, alpha));
}

std::vector<int64_t> continuedFraction(const uint64_t number, const std::size_t k)
{
    std::vector<int64_t> values;
    std::vector<int64_t> alpha;

    double squareNumber = std::sqrt(number);

    uint64_t v_i = 1;
    uint64_t alpha_i = squareNumber;
    uint64_t a_i = std::floor(alpha_i);
    uint64_t u_i = a_i;
    
    for(std::size_t i = 0; i < k + 1; ++i)
    {
        values.push_back(a_i);
        alpha.push_back(alpha_i);

        v_i = (number - (u_i * u_i)) / v_i;
        alpha_i = (squareNumber + u_i) / v_i;
        a_i = std::floor(alpha_i);
        u_i = a_i * v_i - u_i;
    }
    values.push_back(a_i);
    alpha.push_back(alpha_i);

    return values;
}


std::pair<std::vector<int64_t>, bool> factorizeSmoothValue(const int64_t smoothValue, const int64_t number, const std::vector<int64_t>& factorBase)
{
    std::vector<int64_t> exponentsVector(factorBase.size(), 0);

    int64_t smoothValueCopy = smoothValue;
    for(std::size_t i = 1; i < factorBase.size(); ++i)
    {
        while(smoothValueCopy % factorBase[i] == 0)
        {
            smoothValueCopy /= factorBase[i];
            exponentsVector[i] += 1;
        }
    }

    if(smoothValueCopy != 1)
    {
        exponentsVector = std::vector<int64_t>(factorBase.size(), 0);
        smoothValueCopy = number - smoothValue;
        exponentsVector[0] = 1;
        for(std::size_t i = 1; i < factorBase.size(); ++i)
        {
            while(smoothValueCopy % factorBase[i] == 0)
            {
                smoothValueCopy /= factorBase[i];
                exponentsVector[i] += 1;
            }
        }

        if(smoothValueCopy != 1)
        {
            return std::make_pair(std::vector<int64_t>(), true);
        }
    }

    for(auto& it : exponentsVector)
    {
        it %= 2;
    }

    return std::make_pair(exponentsVector, false);
}

std::vector<int64_t> findSmoothValues(const uint64_t number, std::vector<int64_t> factorBase, const std::size_t k)
{
    std::vector<int64_t> smoothValues;

    std::cout << "k: " << k << std::endl;
    auto continuedFractionValues = continuedFraction(number, k);
    std::cout << "continued fraction" << std::endl;

    int64_t b_i = 1;
    int64_t b_ii = 0;

    std::size_t i = 0;
    while(i < k + 1)
    {
        int64_t result = continuedFractionValues[i] * b_i + b_ii;
        int64_t smoothValueCandidate = (result * result) % number;

        auto[factorization, isError] = factorizeSmoothValue(smoothValueCandidate, number, factorBase);

        if(!isError)
        {
            smoothValues.push_back(smoothValueCandidate);
            ++i;
        }

        b_ii = b_i;
        b_i = result;
    }

    return smoothValues;
}

void printSolution(const std::vector<int64_t>& solution) 
{
    std::cout << "solution" << std::endl;
    for(int sol : solution)
        std::cout << sol << " ";
    std::cout << "\n";
}

void xorRows(std::vector<int64_t>& a, std::vector<int64_t>& b)
{
    for(std::size_t i = 0; i < a.size(); ++i)
    {
        b[i] ^= a[i];
    }
}

std::vector<int64_t> gaussGF2(std::vector<std::vector<int64_t>>& A) 
{
    auto matrix = A;
    std::size_t numRows = matrix.size();
    std::size_t numCols = matrix[0].size();
    std::size_t row = 0;
    std::size_t col = 0;

    while (row < numRows && col < numCols) 
    {
        std::size_t pivotRow = row;
        while (pivotRow < numRows && matrix[pivotRow][col] == 0) 
        {
            ++pivotRow;
        }

        if (pivotRow == numRows) 
        {
            ++col;
        } 
        else 
        {
            std::swap(matrix[row], matrix[pivotRow]);

            for (std::size_t i = row + 1; i < numRows; ++i) 
            {
                if (matrix[i][col] == 1) 
                {
                    xorRows(matrix[i], matrix[row]);
                }
            }

            ++row;
            ++col;
        }
    }

    std::vector<int64_t> solution(numCols, 0);
    for (int i = numRows - 1; i >= 0; --i) 
    {
        std::size_t pivotCol = 0;
        while (pivotCol < numCols && matrix[i][pivotCol] == 0) 
        {
            ++pivotCol;
        }
        if (pivotCol < numCols) 
        {
            solution[pivotCol] = matrix[i].back();
            for (std::size_t j = pivotCol + 1; j < numCols - 1; ++j) 
            {
                if (matrix[i][j] == 1) 
                {
                    solution[pivotCol] ^= solution[j];
                }
            }
        }
    }

    return solution;
}

int64_t methodBrilhartMorrison(const uint64_t number)
{
    double alpha = 1.0 / std::sqrt(2);

    while(true)
    {
        auto factorBase = generateFactorBase(number, alpha);
        std::cout << "factor base" << std::endl;

        auto smoothValues = findSmoothValues(number, factorBase, factorBase.size() - 1);
        std::cout << "smooth values" << std::endl;

        std::vector<std::vector<int64_t>> sle;
        for(auto it : smoothValues)
        {
            sle.push_back(factorizeSmoothValue(it, number, factorBase).first);
        }

        std::cout << std::endl;
        auto solution = gaussGF2(sle);
        
        int64_t X = 1;
        for(std::size_t i = 0; i < smoothValues.size(); ++i)
        {
            X *= std::pow(smoothValues[i], solution[i]);
        }

        int64_t Y = 1;
        for(std::size_t j = 0; j < factorBase.size(); ++j)
        {
            int64_t tempResult = 0;
            for(std::size_t i = 0; i < sle.size(); ++i)
            {
                tempResult += solution[i] * sle[i][j];
            }
            tempResult /= 2;

            Y *= std::pow(factorBase[j], tempResult);
        }

        int64_t firstPossibleResult = std::gcd(X + Y, number);
        int64_t secondPossibleResult = std::gcd(X - Y, number);

        if(firstPossibleResult > 1 && firstPossibleResult < number)
        {
            return firstPossibleResult;
        }
        else if(secondPossibleResult > 1 && secondPossibleResult < number)
        {
            return secondPossibleResult;
        }

        alpha += 0.5;
    }
}

void logDivisor(uint64_t divisor, const std::string& method) {
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Divisor found: " << divisor << ", using method: " << method
              << ", at time: " << std::ctime(&now);
}

std::vector<uint64_t> findCanonicalExpansion(uint64_t n) {
    std::vector<uint64_t> factors;
    if (solovey_strassen_test(n)) {
        factors.push_back(n);
        logDivisor(n, "Solovay-Strassen (primality test)");
        return factors;
    }

    auto startTime = std::chrono::system_clock::now();
    while (n > 1) {
        if (solovey_strassen_test(n)) {
            factors.push_back(n);
            logDivisor(n, "Solovay-Strassen (primality test)");
            break;
        }

        uint64_t divisor = methodOfTrialDivisions(n);
        if (divisor > 1 && divisor <= 47) {
            factors.push_back(divisor);
            logDivisor(divisor, "Trial Division");
            n /= divisor;
            continue;
        }

        divisor = rhoPollard(n);
        if (divisor > 1) {
            factors.push_back(divisor);
            logDivisor(divisor, "Pollard's Rho");
            n /= divisor;
            continue;
        }

        divisor = methodBrilhartMorrison(n);
        if (divisor > 1) {
            factors.push_back(divisor);
            logDivisor(divisor, "Brilhart-Morrison");
            n /= divisor;
            continue;
        }

        std::cout << "I can't find the canonical expansion of the number :(" << std::endl;
        break;
    }
    auto endTime = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = endTime - startTime;
    std::cout << "Algorithm execution time: " << elapsedSeconds.count() << "s\n";

    return factors;
}