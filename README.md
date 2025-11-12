# Probability Distribution Calculator

A high-performance C++ probability distribution calculator with interactive parameter input and comprehensive documentation for statistical distributions.

## Features

- **19 Probability Distributions**: Includes both discrete and continuous distributions
- **Special Functions**: Combination (nCr), Permutation (nPr), and Factorial calculations
- **Interactive Interface**: User-friendly menu with detailed parameter descriptions
- **Cumulative Options**: All distributions support CDF (cumulative distribution function)
- **Low-Latency Performance**: Optimized with template metaprogramming and fast mathematical algorithms
- **Built-in Documentation**: Each distribution includes usage examples and when to use them

## Distributions Included

### Discrete Distributions (5)
1. **Binomial Distribution** - Number of successes in fixed trials
2. **Negative Binomial Distribution** - Number of failures before r successes
3. **Poisson Distribution** - Number of events in a fixed interval
4. **Hypergeometric Distribution** - Sampling without replacement
5. **Geometric Distribution** - Trials until first success

### Continuous Distributions (14)
6. **Normal (Gaussian) Distribution** - Bell curve for naturally occurring data
7. **Uniform Distribution** - Equal probability across an interval
8. **Student's t-Distribution** - Small sample hypothesis testing
9. **F-Distribution** - Comparing variances and ANOVA
10. **Gamma Distribution** - Waiting times and continuous positive variables
11. **Chi-Square Distribution** - Goodness-of-fit tests
12. **Exponential Distribution** - Time between events
13. **Beta Distribution** - Probabilities and proportions (0 to 1)
14. **Weibull Distribution** - Reliability analysis and failure rates
15. **Lognormal Distribution** - Stock prices and multiplicative processes
16. **Cauchy Distribution** - Heavy-tailed phenomena and physics
17. **Pareto Distribution** - Wealth distribution and 80-20 rule
18. **Laplace Distribution** - Financial returns and L1 regularization
19. **Logistic Distribution** - Growth models and logistic regression

### Special Functions (3)
20. **Combination (nCr)** - Choose r items from n (order doesn't matter)
21. **Permutation (nPr)** - Arrange r items from n (order matters)
22. **Factorial (n!)** - Product of all positive integers up to n

## Building the Project

### Prerequisites
- C++17 compatible compiler (g++, clang++)
- CMake 3.10 or higher

### Build Instructions

```bash
# Create build directory
mkdir build
cd build

# Configure and build
cmake ..
make

# Run the program
./prob_calc
```

## Usage Examples

### Example 1: Binomial Distribution (Coin Flipping)
```
Question: If you flip a fair coin 10 times, what's the probability of getting exactly 6 heads?

Input:
- n = 10 (number of trials)
- p = 0.5 (probability of success)
- k = 6 (number of successes)

Output:
- P(X = 6) = 0.205078 (exactly 6 heads)
```

### Example 2: Normal Distribution (Heights)
```
Question: If average height is 170 cm with standard deviation 10 cm,
what's the probability someone is between 165 and 175 cm?

Input:
- μ = 170 (mean)
- σ = 10 (standard deviation)
- Calculate area between 165 and 175

Output:
- P(165 ≤ X ≤ 175) = 0.383 (38.3% of people)
```

### Example 3: Combination (Lottery)
```
Question: How many ways can you choose 6 numbers from 49?

Input:
- n = 49
- r = 6

Output:
- C(49, 6) = 13,983,816 combinations
```

## Performance Optimizations

This calculator uses several optimization techniques:

1. **Template Metaprogramming**: Compile-time factorial calculations for small values
2. **Memoization**: Caching factorial results to avoid recalculation
3. **Logarithmic Calculations**: Using log-space for numerical stability
4. **Fast Mathematical Functions**: Optimized implementations of gamma, beta, and error functions
5. **Compiler Optimizations**: -O3, -march=native, -ffast-math flags

## Mathematical Algorithms

- **Gamma Function**: Lanczos approximation for high accuracy
- **Beta Function**: Computed using gamma function ratios
- **Incomplete Beta**: Continued fractions method
- **Error Function**: Polynomial approximation
- **Normal CDF**: Computed using error function

## Project Structure

```
probaillty_calculator/
├── main.cpp                        # Main program with interactive menu
├── math_utils.hpp                  # Core mathematical utilities
├── discrete_distributions.hpp      # Discrete probability distributions
├── continuous_distributions.hpp    # Continuous probability distributions
├── special_functions.hpp           # nCr, nPr, factorial
├── CMakeLists.txt                  # Build configuration
└── README.md                       # This file
```

## When to Use Each Distribution

### Discrete Distributions

- **Binomial**: Fixed number of trials, two outcomes (success/failure), independent trials
- **Negative Binomial**: Counting failures before achieving fixed successes
- **Poisson**: Counting events over time/space when events are rare
- **Hypergeometric**: Sampling without replacement from finite population
- **Geometric**: Waiting for first success

### Continuous Distributions

- **Normal**: Natural phenomena, measurement errors, Central Limit Theorem
- **Uniform**: Random number generation, equal likelihood
- **Student's t**: Small samples, unknown population variance
- **F-Distribution**: Comparing two variances, ANOVA
- **Gamma**: Waiting times for multiple events, modeling positive continuous data
- **Chi-Square**: Goodness-of-fit tests, variance testing
- **Exponential**: Time between events in Poisson process
- **Beta**: Modeling probabilities and proportions
- **Weibull**: Reliability analysis, failure rates, wind speed modeling
- **Lognormal**: Stock prices, income distributions, multiplicative processes
- **Cauchy**: Heavy-tailed phenomena, ratio of normal variables, physics
- **Pareto**: Wealth distribution, 80-20 rule, insurance claims
- **Laplace**: Financial returns (fat tails), signal processing, L1 regularization
- **Logistic**: Growth models, technology adoption, logistic regression

## Contributing

Feel free to submit issues or pull requests for:
- Additional distributions
- Performance improvements
- Bug fixes
- Documentation enhancements

## License

This project is open source and available for educational and commercial use.

## Author

Created as a comprehensive probability calculator for statistical analysis and education.
