#pragma once

#include <cmath>
#include <limits>
#include <algorithm>

namespace ProbCalc {

// Compile-time factorial using template metaprogramming
template<int N>
struct Factorial {
    static constexpr long long value = N * Factorial<N - 1>::value;
};

template<>
struct Factorial<0> {
    static constexpr long long value = 1;
};

// Runtime factorial for large numbers (with memoization)
inline double factorial(int n) {
    static double cache[171] = {0}; // double can hold up to 170!
    if (n < 0) return std::numeric_limits<double>::quiet_NaN();
    if (n > 170) return std::numeric_limits<double>::infinity();

    if (cache[n] > 0) return cache[n];

    if (n == 0 || n == 1) {
        cache[n] = 1.0;
        return 1.0;
    }

    cache[n] = n * factorial(n - 1);
    return cache[n];
}

// Log factorial for better numerical stability
inline double log_factorial(int n) {
    if (n < 0) return std::numeric_limits<double>::quiet_NaN();
    if (n == 0 || n == 1) return 0.0;

    double result = 0.0;
    for (int i = 2; i <= n; ++i) {
        result += std::log(static_cast<double>(i));
    }
    return result;
}

// Combination nCr using template metaprogramming where possible
template<int N, int R>
struct Combination {
    static constexpr long long value =
        Factorial<N>::value / (Factorial<R>::value * Factorial<N - R>::value);
};

// Runtime combination
inline double combination(int n, int r) {
    if (r < 0 || r > n) return 0.0;
    if (r == 0 || r == n) return 1.0;

    // Use symmetry: C(n,r) = C(n, n-r)
    if (r > n - r) r = n - r;

    // Use log to avoid overflow
    return std::exp(log_factorial(n) - log_factorial(r) - log_factorial(n - r));
}

// Permutation nPr
inline double permutation(int n, int r) {
    if (r < 0 || r > n) return 0.0;
    if (r == 0) return 1.0;

    return std::exp(log_factorial(n) - log_factorial(n - r));
}

// Gamma function using Lanczos approximation
inline double gamma(double z) {
    if (z < 0.5) {
        // Reflection formula: Gamma(z) * Gamma(1-z) = pi / sin(pi*z)
        return M_PI / (std::sin(M_PI * z) * gamma(1.0 - z));
    }

    // Lanczos coefficients
    static const double g = 7.0;
    static const double coef[] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    z -= 1.0;
    double x = coef[0];
    for (int i = 1; i < 9; ++i) {
        x += coef[i] / (z + i);
    }

    double t = z + g + 0.5;
    return std::sqrt(2.0 * M_PI) * std::pow(t, z + 0.5) * std::exp(-t) * x;
}

// Log gamma function
inline double log_gamma(double z) {
    if (z <= 0) return std::numeric_limits<double>::quiet_NaN();
    return std::log(gamma(z));
}

// Beta function: B(a,b) = Gamma(a) * Gamma(b) / Gamma(a+b)
inline double beta(double a, double b) {
    return std::exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b));
}

// Error function using polynomial approximation
inline double erf(double x) {
    // Constants for approximation
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    int sign = (x < 0) ? -1 : 1;
    x = std::abs(x);

    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

    return sign * y;
}

// Complementary error function
inline double erfc(double x) {
    return 1.0 - erf(x);
}

// Incomplete beta function using continued fractions
inline double incomplete_beta(double x, double a, double b, int max_iter = 200) {
    if (x < 0.0 || x > 1.0) return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return 0.0;
    if (x == 1.0) return 1.0;

    // Use symmetry relation if needed
    bool flag = false;
    if (x > (a + 1.0) / (a + b + 2.0)) {
        flag = true;
        std::swap(a, b);
        x = 1.0 - x;
    }

    double lbeta = log_gamma(a + b) - log_gamma(a) - log_gamma(b) +
                   a * std::log(x) + b * std::log(1.0 - x);

    // Continued fraction
    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;

    if (std::abs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= max_iter; ++m) {
        int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (std::abs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (std::abs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double del = d * c;
        h *= del;

        if (std::abs(del - 1.0) < 1e-10) break;
    }

    double result = std::exp(lbeta) * h / a;
    return flag ? 1.0 - result : result;
}

// Normal CDF using error function
inline double normal_cdf(double x, double mu = 0.0, double sigma = 1.0) {
    return 0.5 * (1.0 + erf((x - mu) / (sigma * std::sqrt(2.0))));
}

// Normal PDF
inline double normal_pdf(double x, double mu = 0.0, double sigma = 1.0) {
    double z = (x - mu) / sigma;
    return (1.0 / (sigma * std::sqrt(2.0 * M_PI))) * std::exp(-0.5 * z * z);
}

} // namespace ProbCalc
