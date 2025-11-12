#pragma once

#include "math_utils.hpp"
#include <cmath>
#include <algorithm>
#include <string>

namespace ProbCalc {

// Binomial Distribution
class BinomialDistribution {
public:
    static std::string get_name() { return "Binomial Distribution"; }

    static std::string get_description() {
        return "Used for: Number of successes in n independent trials with probability p.\n"
               "Example: Flipping a coin n times, counting heads.\n"
               "Example: Number of defective items in a batch.\n"
               "Example: Number of successful sales calls out of n attempts.";
    }

    static double pmf(int k, int n, double p) {
        if (k < 0 || k > n || p < 0.0 || p > 1.0) return 0.0;
        if (p == 0.0) return (k == 0) ? 1.0 : 0.0;
        if (p == 1.0) return (k == n) ? 1.0 : 0.0;

        // Use log for numerical stability
        double log_prob = log_factorial(n) - log_factorial(k) - log_factorial(n - k) +
                         k * std::log(p) + (n - k) * std::log(1.0 - p);
        return std::exp(log_prob);
    }

    static double cdf(int k, int n, double p) {
        if (k < 0) return 0.0;
        if (k >= n) return 1.0;

        double sum = 0.0;
        for (int i = 0; i <= k; ++i) {
            sum += pmf(i, n, p);
        }
        return sum;
    }
};

// Negative Binomial Distribution
class NegativeBinomialDistribution {
public:
    static std::string get_name() { return "Negative Binomial Distribution"; }

    static std::string get_description() {
        return "Used for: Number of failures before r successes occur.\n"
               "Example: Number of coin flips before getting 5 heads.\n"
               "Example: Number of items inspected before finding 3 defects.\n"
               "Example: Over-dispersed count data (alternative to Poisson).";
    }

    static double pmf(int k, int r, double p) {
        if (k < 0 || r <= 0 || p <= 0.0 || p > 1.0) return 0.0;

        // P(X = k) = C(k+r-1, k) * p^r * (1-p)^k
        double log_prob = log_factorial(k + r - 1) - log_factorial(k) - log_factorial(r - 1) +
                         r * std::log(p) + k * std::log(1.0 - p);
        return std::exp(log_prob);
    }

    static double cdf(int k, int r, double p) {
        if (k < 0) return 0.0;

        double sum = 0.0;
        for (int i = 0; i <= k; ++i) {
            sum += pmf(i, r, p);
        }
        return sum;
    }
};

// Poisson Distribution
class PoissonDistribution {
public:
    static std::string get_name() { return "Poisson Distribution"; }

    static std::string get_description() {
        return "Used for: Number of events occurring in a fixed interval.\n"
               "Example: Number of customers arriving at a store per hour.\n"
               "Example: Number of phone calls received in a day.\n"
               "Example: Number of defects in a unit area of material.";
    }

    static double pmf(int k, double lambda) {
        if (k < 0 || lambda <= 0.0) return 0.0;

        // P(X = k) = (lambda^k * e^(-lambda)) / k!
        double log_prob = k * std::log(lambda) - lambda - log_factorial(k);
        return std::exp(log_prob);
    }

    static double cdf(int k, double lambda) {
        if (k < 0) return 0.0;

        double sum = 0.0;
        for (int i = 0; i <= k; ++i) {
            sum += pmf(i, lambda);
        }
        return sum;
    }
};

// Hypergeometric Distribution
class HypergeometricDistribution {
public:
    static std::string get_name() { return "Hypergeometric Distribution"; }

    static std::string get_description() {
        return "Used for: Sampling without replacement from a finite population.\n"
               "Example: Drawing cards from a deck without replacement.\n"
               "Example: Quality control - selecting items from a batch.\n"
               "Example: Number of defective items in a sample from a lot.";
    }

    static double pmf(int k, int N, int K, int n) {
        // N = population size
        // K = number of success states in population
        // n = number of draws
        // k = number of observed successes

        if (k < 0 || k > n || k > K || (n - k) > (N - K)) return 0.0;

        // P(X = k) = C(K, k) * C(N-K, n-k) / C(N, n)
        double log_prob = (log_factorial(K) - log_factorial(k) - log_factorial(K - k)) +
                         (log_factorial(N - K) - log_factorial(n - k) - log_factorial(N - K - n + k)) -
                         (log_factorial(N) - log_factorial(n) - log_factorial(N - n));
        return std::exp(log_prob);
    }

    static double cdf(int k, int N, int K, int n) {
        if (k < 0) return 0.0;

        double sum = 0.0;
        int max_k = std::min(n, K);
        for (int i = 0; i <= std::min(k, max_k); ++i) {
            sum += pmf(i, N, K, n);
        }
        return sum;
    }
};

// Geometric Distribution (special case of negative binomial with r=1)
class GeometricDistribution {
public:
    static std::string get_name() { return "Geometric Distribution"; }

    static std::string get_description() {
        return "Used for: Number of trials until the first success.\n"
               "Example: Number of coin flips until first heads.\n"
               "Example: Number of products tested until finding first defect.\n"
               "Example: Number of attempts until first successful connection.";
    }

    static double pmf(int k, double p) {
        if (k < 1 || p <= 0.0 || p > 1.0) return 0.0;

        // P(X = k) = (1-p)^(k-1) * p
        return std::pow(1.0 - p, k - 1) * p;
    }

    static double cdf(int k, double p) {
        if (k < 1) return 0.0;

        // P(X <= k) = 1 - (1-p)^k
        return 1.0 - std::pow(1.0 - p, k);
    }
};

} // namespace ProbCalc
