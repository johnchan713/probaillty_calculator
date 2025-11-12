#pragma once

#include "math_utils.hpp"
#include <cmath>
#include <algorithm>
#include <string>

namespace ProbCalc {

// Normal (Gaussian) Distribution
class NormalDistribution {
public:
    static std::string get_name() { return "Normal (Gaussian) Distribution"; }

    static std::string get_description() {
        return "Used for: Continuous data that clusters around a mean.\n"
               "Example: Heights, weights, test scores, measurement errors.\n"
               "Example: Natural phenomena (IQ scores, blood pressure).\n"
               "Example: Central limit theorem - sums of random variables.";
    }

    static double pdf(double x, double mu, double sigma) {
        if (sigma <= 0.0) return 0.0;
        return normal_pdf(x, mu, sigma);
    }

    static double cdf(double x, double mu, double sigma) {
        if (sigma <= 0.0) return 0.0;
        return normal_cdf(x, mu, sigma);
    }

    static double area_between(double a, double b, double mu, double sigma) {
        return cdf(b, mu, sigma) - cdf(a, mu, sigma);
    }
};

// Uniform Distribution
class UniformDistribution {
public:
    static std::string get_name() { return "Uniform Distribution"; }

    static std::string get_description() {
        return "Used for: All outcomes equally likely in an interval.\n"
               "Example: Random number generation between a and b.\n"
               "Example: Round-off errors in measurements.\n"
               "Example: Arrival time within a given time window.";
    }

    static double pdf(double x, double a, double b) {
        if (a >= b) return 0.0;
        if (x < a || x > b) return 0.0;
        return 1.0 / (b - a);
    }

    static double cdf(double x, double a, double b) {
        if (a >= b) return 0.0;
        if (x < a) return 0.0;
        if (x > b) return 1.0;
        return (x - a) / (b - a);
    }

    static double area_between(double x1, double x2, double a, double b) {
        return cdf(x2, a, b) - cdf(x1, a, b);
    }
};

// Student's t-Distribution
class StudentTDistribution {
public:
    static std::string get_name() { return "Student's t-Distribution"; }

    static std::string get_description() {
        return "Used for: Small sample statistics and hypothesis testing.\n"
               "Example: Confidence intervals with unknown population variance.\n"
               "Example: T-tests for comparing means.\n"
               "Example: Regression analysis with small samples.";
    }

    static double pdf(double x, double df) {
        if (df <= 0.0) return 0.0;

        double num = gamma((df + 1.0) / 2.0);
        double den = std::sqrt(df * M_PI) * gamma(df / 2.0);
        double factor = std::pow(1.0 + (x * x) / df, -(df + 1.0) / 2.0);

        return (num / den) * factor;
    }

    static double cdf(double x, double df) {
        if (df <= 0.0) return 0.0;

        double t2 = x * x / df;
        double ival = incomplete_beta(df / (df + t2), df / 2.0, 0.5);

        if (x >= 0.0) {
            return 0.5 + 0.5 * (1.0 - ival);
        } else {
            return 0.5 * ival;
        }
    }

    static double area_between(double a, double b, double df) {
        return cdf(b, df) - cdf(a, df);
    }
};

// F-Distribution
class FDistribution {
public:
    static std::string get_name() { return "F-Distribution"; }

    static std::string get_description() {
        return "Used for: Comparing variances and ANOVA.\n"
               "Example: Testing equality of variances from two populations.\n"
               "Example: Analysis of variance (ANOVA) tests.\n"
               "Example: Regression analysis (F-test for overall significance).";
    }

    static double pdf(double x, double d1, double d2) {
        if (x <= 0.0 || d1 <= 0.0 || d2 <= 0.0) return 0.0;

        double num = std::sqrt(std::pow(d1 * x, d1) * std::pow(d2, d2) /
                              std::pow(d1 * x + d2, d1 + d2));
        double den = x * beta(d1 / 2.0, d2 / 2.0);

        return num / den;
    }

    static double cdf(double x, double d1, double d2) {
        if (x <= 0.0) return 0.0;
        if (d1 <= 0.0 || d2 <= 0.0) return 0.0;

        double ival = d1 * x / (d1 * x + d2);
        return incomplete_beta(ival, d1 / 2.0, d2 / 2.0);
    }

    static double area_between(double a, double b, double d1, double d2) {
        return cdf(b, d1, d2) - cdf(a, d1, d2);
    }
};

// Gamma Distribution
class GammaDistribution {
public:
    static std::string get_name() { return "Gamma Distribution"; }

    static std::string get_description() {
        return "Used for: Waiting times and continuous positive variables.\n"
               "Example: Time until k events occur (k = shape parameter).\n"
               "Example: Rainfall amounts, insurance claims.\n"
               "Example: Time between events in a Poisson process.";
    }

    static double pdf(double x, double shape, double scale) {
        if (x <= 0.0 || shape <= 0.0 || scale <= 0.0) return 0.0;

        double log_pdf = (shape - 1.0) * std::log(x) - x / scale -
                        shape * std::log(scale) - log_gamma(shape);
        return std::exp(log_pdf);
    }

    static double cdf(double x, double shape, double scale) {
        if (x <= 0.0) return 0.0;
        if (shape <= 0.0 || scale <= 0.0) return 0.0;

        // Use incomplete gamma function approximation
        double z = x / scale;
        return incomplete_beta(z / (z + shape), shape, 1.0);
    }

    static double area_between(double a, double b, double shape, double scale) {
        return cdf(b, shape, scale) - cdf(a, shape, scale);
    }
};

// Chi-Square Distribution
class ChiSquareDistribution {
public:
    static std::string get_name() { return "Chi-Square Distribution"; }

    static std::string get_description() {
        return "Used for: Goodness-of-fit tests and variance estimation.\n"
               "Example: Chi-square test for independence in contingency tables.\n"
               "Example: Goodness-of-fit tests for distributions.\n"
               "Example: Confidence intervals for variance.\n"
               "Note: Chi-square is a special case of Gamma with shape=k/2, scale=2.";
    }

    static double pdf(double x, double k) {
        if (x <= 0.0 || k <= 0.0) return 0.0;

        // Chi-square is Gamma(k/2, 2)
        double shape = k / 2.0;
        double log_pdf = (shape - 1.0) * std::log(x) - x / 2.0 -
                        shape * std::log(2.0) - log_gamma(shape);
        return std::exp(log_pdf);
    }

    static double cdf(double x, double k) {
        if (x <= 0.0) return 0.0;
        if (k <= 0.0) return 0.0;

        // Use incomplete gamma
        double z = x / 2.0;
        double shape = k / 2.0;
        return incomplete_beta(z / (z + shape), shape, 1.0);
    }

    static double area_between(double a, double b, double k) {
        return cdf(b, k) - cdf(a, k);
    }
};

// Exponential Distribution
class ExponentialDistribution {
public:
    static std::string get_name() { return "Exponential Distribution"; }

    static std::string get_description() {
        return "Used for: Time between events in a Poisson process.\n"
               "Example: Time between arrivals of customers.\n"
               "Example: Lifetime of electronic components.\n"
               "Example: Time until next earthquake or accident.";
    }

    static double pdf(double x, double lambda) {
        if (x < 0.0 || lambda <= 0.0) return 0.0;
        return lambda * std::exp(-lambda * x);
    }

    static double cdf(double x, double lambda) {
        if (x < 0.0) return 0.0;
        if (lambda <= 0.0) return 0.0;
        return 1.0 - std::exp(-lambda * x);
    }

    static double area_between(double a, double b, double lambda) {
        return cdf(b, lambda) - cdf(a, lambda);
    }
};

// Beta Distribution
class BetaDistribution {
public:
    static std::string get_name() { return "Beta Distribution"; }

    static std::string get_description() {
        return "Used for: Probabilities and proportions (bounded between 0 and 1).\n"
               "Example: Modeling probabilities of success rates.\n"
               "Example: Bayesian prior for binomial probability.\n"
               "Example: Project completion percentages.";
    }

    static double pdf(double x, double alpha, double beta_param) {
        if (x <= 0.0 || x >= 1.0 || alpha <= 0.0 || beta_param <= 0.0) return 0.0;

        double log_pdf = (alpha - 1.0) * std::log(x) + (beta_param - 1.0) * std::log(1.0 - x) -
                        log_gamma(alpha) - log_gamma(beta_param) + log_gamma(alpha + beta_param);
        return std::exp(log_pdf);
    }

    static double cdf(double x, double alpha, double beta_param) {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;
        if (alpha <= 0.0 || beta_param <= 0.0) return 0.0;

        return incomplete_beta(x, alpha, beta_param);
    }

    static double area_between(double a, double b, double alpha, double beta_param) {
        return cdf(b, alpha, beta_param) - cdf(a, alpha, beta_param);
    }
};

// Weibull Distribution
class WeibullDistribution {
public:
    static std::string get_name() { return "Weibull Distribution"; }

    static std::string get_description() {
        return "Used for: Reliability analysis, failure rates, and lifetime modeling.\n"
               "Example: Time-to-failure of mechanical components.\n"
               "Example: Wind speed distributions in renewable energy.\n"
               "Example: Product lifetime and warranty analysis.\n"
               "Parameters: k (shape, k>1 increasing hazard), λ (scale, characteristic life).";
    }

    static double pdf(double x, double k, double lambda) {
        if (x < 0.0 || k <= 0.0 || lambda <= 0.0) return 0.0;

        double xk = std::pow(x / lambda, k);
        return (k / lambda) * std::pow(x / lambda, k - 1.0) * std::exp(-xk);
    }

    static double cdf(double x, double k, double lambda) {
        if (x < 0.0) return 0.0;
        if (k <= 0.0 || lambda <= 0.0) return 0.0;

        return 1.0 - std::exp(-std::pow(x / lambda, k));
    }

    static double area_between(double a, double b, double k, double lambda) {
        return cdf(b, k, lambda) - cdf(a, k, lambda);
    }
};

// Lognormal Distribution
class LognormalDistribution {
public:
    static std::string get_name() { return "Lognormal Distribution"; }

    static std::string get_description() {
        return "Used for: Positive variables that are products of many factors.\n"
               "Example: Stock prices, asset returns, income distributions.\n"
               "Example: Size of biological organisms, particle sizes.\n"
               "Example: Time to complete tasks, file sizes.\n"
               "Parameters: μ (mean of log(X)), σ (std dev of log(X)).";
    }

    static double pdf(double x, double mu, double sigma) {
        if (x <= 0.0 || sigma <= 0.0) return 0.0;

        double log_x = std::log(x);
        double z = (log_x - mu) / sigma;
        return (1.0 / (x * sigma * std::sqrt(2.0 * M_PI))) *
               std::exp(-0.5 * z * z);
    }

    static double cdf(double x, double mu, double sigma) {
        if (x <= 0.0) return 0.0;
        if (sigma <= 0.0) return 0.0;

        double log_x = std::log(x);
        return normal_cdf(log_x, mu, sigma);
    }

    static double area_between(double a, double b, double mu, double sigma) {
        return cdf(b, mu, sigma) - cdf(a, mu, sigma);
    }
};

// Cauchy Distribution
class CauchyDistribution {
public:
    static std::string get_name() { return "Cauchy Distribution"; }

    static std::string get_description() {
        return "Used for: Heavy-tailed phenomena, resonance in physics.\n"
               "Example: Ratio of two normal random variables.\n"
               "Example: Spectral line shapes in physics.\n"
               "Example: Modeling outliers and extreme values.\n"
               "Parameters: x₀ (location/median), γ (scale/half-width at half-maximum).";
    }

    static double pdf(double x, double x0, double gamma) {
        if (gamma <= 0.0) return 0.0;

        double z = (x - x0) / gamma;
        return 1.0 / (M_PI * gamma * (1.0 + z * z));
    }

    static double cdf(double x, double x0, double gamma) {
        if (gamma <= 0.0) return 0.0;

        return 0.5 + (1.0 / M_PI) * std::atan((x - x0) / gamma);
    }

    static double area_between(double a, double b, double x0, double gamma) {
        return cdf(b, x0, gamma) - cdf(a, x0, gamma);
    }
};

// Pareto Distribution
class ParetoDistribution {
public:
    static std::string get_name() { return "Pareto Distribution"; }

    static std::string get_description() {
        return "Used for: Wealth distribution, 80-20 rule phenomena.\n"
               "Example: Income distribution (top earners).\n"
               "Example: City population sizes.\n"
               "Example: Insurance claims, large losses.\n"
               "Parameters: xₘ (minimum value), α (shape, controls tail heaviness).";
    }

    static double pdf(double x, double xm, double alpha) {
        if (x < xm || xm <= 0.0 || alpha <= 0.0) return 0.0;

        return (alpha * std::pow(xm, alpha)) / std::pow(x, alpha + 1.0);
    }

    static double cdf(double x, double xm, double alpha) {
        if (x < xm) return 0.0;
        if (xm <= 0.0 || alpha <= 0.0) return 0.0;

        return 1.0 - std::pow(xm / x, alpha);
    }

    static double area_between(double a, double b, double xm, double alpha) {
        return cdf(b, xm, alpha) - cdf(a, xm, alpha);
    }
};

// Laplace Distribution
class LaplaceDistribution {
public:
    static std::string get_name() { return "Laplace Distribution"; }

    static std::string get_description() {
        return "Used for: Modeling differences, L1 regularization, signal processing.\n"
               "Example: Financial returns (fatter tails than normal).\n"
               "Example: Difference between two exponential variables.\n"
               "Example: Lasso regression (L1 penalty).\n"
               "Parameters: μ (location/mean/median), b (scale/diversity).";
    }

    static double pdf(double x, double mu, double b) {
        if (b <= 0.0) return 0.0;

        return (1.0 / (2.0 * b)) * std::exp(-std::abs(x - mu) / b);
    }

    static double cdf(double x, double mu, double b) {
        if (b <= 0.0) return 0.0;

        if (x < mu) {
            return 0.5 * std::exp((x - mu) / b);
        } else {
            return 1.0 - 0.5 * std::exp(-(x - mu) / b);
        }
    }

    static double area_between(double a, double b_upper, double mu, double b) {
        return cdf(b_upper, mu, b) - cdf(a, mu, b);
    }
};

// Logistic Distribution
class LogisticDistribution {
public:
    static std::string get_name() { return "Logistic Distribution"; }

    static std::string get_description() {
        return "Used for: Growth models, logistic regression, neural networks.\n"
               "Example: Population growth (S-curve).\n"
               "Example: Adoption of new technology.\n"
               "Example: Logistic regression classification.\n"
               "Parameters: μ (location/mean), s (scale, proportional to std dev).";
    }

    static double pdf(double x, double mu, double s) {
        if (s <= 0.0) return 0.0;

        double z = (x - mu) / s;
        double exp_z = std::exp(-z);
        return exp_z / (s * std::pow(1.0 + exp_z, 2.0));
    }

    static double cdf(double x, double mu, double s) {
        if (s <= 0.0) return 0.0;

        double z = (x - mu) / s;
        return 1.0 / (1.0 + std::exp(-z));
    }

    static double area_between(double a, double b, double mu, double s) {
        return cdf(b, mu, s) - cdf(a, mu, s);
    }
};

} // namespace ProbCalc
