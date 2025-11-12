#include "discrete_distributions.hpp"
#include "continuous_distributions.hpp"
#include "special_functions.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

using namespace ProbCalc;

// Helper function to get valid integer input
int get_int_input(const std::string& prompt) {
    int value;
    while (true) {
        std::cout << prompt;
        if (std::cin >> value) {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return value;
        } else {
            std::cout << "Invalid input. Please enter a valid integer.\n";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
}

// Helper function to get valid double input
double get_double_input(const std::string& prompt) {
    double value;
    while (true) {
        std::cout << prompt;
        if (std::cin >> value) {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return value;
        } else {
            std::cout << "Invalid input. Please enter a valid number.\n";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
}

// Helper function to display menu header
void display_header(const std::string& title) {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "  " << title << "\n";
    std::cout << "========================================\n\n";
}

// Function to handle binomial distribution
void handle_binomial() {
    display_header(BinomialDistribution::get_name());
    std::cout << BinomialDistribution::get_description() << "\n\n";

    int n = get_int_input("Enter n (number of trials, e.g., how many times you flip a coin): ");
    double p = get_double_input("Enter p (probability of success, e.g., 0.5 for fair coin): ");
    int k = get_int_input("Enter k (number of successes you want to calculate, e.g., 3 heads): ");

    double pmf = BinomialDistribution::pmf(k, n, p);
    double cdf = BinomialDistribution::cdf(k, n, p);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResults:\n";
    std::cout << "  P(X = " << k << ") = " << pmf << " (exactly " << k << " successes)\n";
    std::cout << "  P(X <= " << k << ") = " << cdf << " (at most " << k << " successes)\n";
    std::cout << "  P(X > " << k << ") = " << (1.0 - cdf) << " (more than " << k << " successes)\n";
}

// Function to handle negative binomial distribution
void handle_negative_binomial() {
    display_header(NegativeBinomialDistribution::get_name());
    std::cout << NegativeBinomialDistribution::get_description() << "\n\n";

    int r = get_int_input("Enter r (number of successes required, e.g., 5 heads): ");
    double p = get_double_input("Enter p (probability of success on each trial, e.g., 0.5): ");
    int k = get_int_input("Enter k (number of failures before r successes): ");

    double pmf = NegativeBinomialDistribution::pmf(k, r, p);
    double cdf = NegativeBinomialDistribution::cdf(k, r, p);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResults:\n";
    std::cout << "  P(X = " << k << ") = " << pmf << " (exactly " << k << " failures)\n";
    std::cout << "  P(X <= " << k << ") = " << cdf << " (at most " << k << " failures)\n";
}

// Function to handle Poisson distribution
void handle_poisson() {
    display_header(PoissonDistribution::get_name());
    std::cout << PoissonDistribution::get_description() << "\n\n";

    double lambda = get_double_input("Enter λ (lambda - average rate, e.g., 5 customers per hour): ");
    int k = get_int_input("Enter k (number of events, e.g., 7 customers): ");

    double pmf = PoissonDistribution::pmf(k, lambda);
    double cdf = PoissonDistribution::cdf(k, lambda);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResults:\n";
    std::cout << "  P(X = " << k << ") = " << pmf << " (exactly " << k << " events)\n";
    std::cout << "  P(X <= " << k << ") = " << cdf << " (at most " << k << " events)\n";
}

// Function to handle hypergeometric distribution
void handle_hypergeometric() {
    display_header(HypergeometricDistribution::get_name());
    std::cout << HypergeometricDistribution::get_description() << "\n\n";

    int N = get_int_input("Enter N (total population size, e.g., 52 cards in deck): ");
    int K = get_int_input("Enter K (number of success states in population, e.g., 13 hearts): ");
    int n = get_int_input("Enter n (number of draws/samples, e.g., 5 cards drawn): ");
    int k = get_int_input("Enter k (number of successes desired, e.g., 2 hearts): ");

    double pmf = HypergeometricDistribution::pmf(k, N, K, n);
    double cdf = HypergeometricDistribution::cdf(k, N, K, n);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResults:\n";
    std::cout << "  P(X = " << k << ") = " << pmf << " (exactly " << k << " successes)\n";
    std::cout << "  P(X <= " << k << ") = " << cdf << " (at most " << k << " successes)\n";
}

// Function to handle geometric distribution
void handle_geometric() {
    display_header(GeometricDistribution::get_name());
    std::cout << GeometricDistribution::get_description() << "\n\n";

    double p = get_double_input("Enter p (probability of success on each trial, e.g., 0.5): ");
    int k = get_int_input("Enter k (number of trials until first success, e.g., 3): ");

    double pmf = GeometricDistribution::pmf(k, p);
    double cdf = GeometricDistribution::cdf(k, p);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nResults:\n";
    std::cout << "  P(X = " << k << ") = " << pmf << " (first success on trial " << k << ")\n";
    std::cout << "  P(X <= " << k << ") = " << cdf << " (first success within " << k << " trials)\n";
}

// Function to handle normal distribution
void handle_normal() {
    display_header(NormalDistribution::get_name());
    std::cout << NormalDistribution::get_description() << "\n\n";

    double mu = get_double_input("Enter μ (mu - mean, e.g., average height 170 cm): ");
    double sigma = get_double_input("Enter σ (sigma - standard deviation, e.g., 10 cm): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x (the value, e.g., 180 cm): ");
        double pdf = NormalDistribution::pdf(x, mu, sigma);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x (upper bound): ");
        double cdf = NormalDistribution::cdf(x, mu, sigma);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double a = get_double_input("Enter lower bound a: ");
        double b = get_double_input("Enter upper bound b: ");
        double area = NormalDistribution::area_between(a, b, mu, sigma);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << a << " <= X <= " << b << ") = " << area << "\n";
    }
}

// Function to handle uniform distribution
void handle_uniform() {
    display_header(UniformDistribution::get_name());
    std::cout << UniformDistribution::get_description() << "\n\n";

    double a = get_double_input("Enter a (lower bound of interval): ");
    double b = get_double_input("Enter b (upper bound of interval): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x: ");
        double pdf = UniformDistribution::pdf(x, a, b);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x: ");
        double cdf = UniformDistribution::cdf(x, a, b);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double x1 = get_double_input("Enter lower bound: ");
        double x2 = get_double_input("Enter upper bound: ");
        double area = UniformDistribution::area_between(x1, x2, a, b);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << x1 << " <= X <= " << x2 << ") = " << area << "\n";
    }
}

// Function to handle Student's t-distribution
void handle_student_t() {
    display_header(StudentTDistribution::get_name());
    std::cout << StudentTDistribution::get_description() << "\n\n";

    double df = get_double_input("Enter degrees of freedom (df, e.g., n-1 for sample size n): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x (t-statistic): ");
        double pdf = StudentTDistribution::pdf(x, df);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x (t-statistic): ");
        double cdf = StudentTDistribution::cdf(x, df);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double a = get_double_input("Enter lower bound: ");
        double b = get_double_input("Enter upper bound: ");
        double area = StudentTDistribution::area_between(a, b, df);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << a << " <= X <= " << b << ") = " << area << "\n";
    }
}

// Function to handle F-distribution
void handle_f_distribution() {
    display_header(FDistribution::get_name());
    std::cout << FDistribution::get_description() << "\n\n";

    double d1 = get_double_input("Enter d1 (degrees of freedom for numerator): ");
    double d2 = get_double_input("Enter d2 (degrees of freedom for denominator): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x (F-statistic): ");
        double pdf = FDistribution::pdf(x, d1, d2);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x (F-statistic): ");
        double cdf = FDistribution::cdf(x, d1, d2);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double a = get_double_input("Enter lower bound: ");
        double b = get_double_input("Enter upper bound: ");
        double area = FDistribution::area_between(a, b, d1, d2);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << a << " <= X <= " << b << ") = " << area << "\n";
    }
}

// Function to handle Gamma distribution
void handle_gamma() {
    display_header(GammaDistribution::get_name());
    std::cout << GammaDistribution::get_description() << "\n\n";

    double shape = get_double_input("Enter shape parameter k (k > 0): ");
    double scale = get_double_input("Enter scale parameter θ (theta > 0): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x: ");
        double pdf = GammaDistribution::pdf(x, shape, scale);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x: ");
        double cdf = GammaDistribution::cdf(x, shape, scale);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double a = get_double_input("Enter lower bound: ");
        double b = get_double_input("Enter upper bound: ");
        double area = GammaDistribution::area_between(a, b, shape, scale);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << a << " <= X <= " << b << ") = " << area << "\n";
    }
}

// Function to handle Chi-square distribution
void handle_chi_square() {
    display_header(ChiSquareDistribution::get_name());
    std::cout << ChiSquareDistribution::get_description() << "\n\n";

    double k = get_double_input("Enter k (degrees of freedom): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x (chi-square statistic): ");
        double pdf = ChiSquareDistribution::pdf(x, k);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x (chi-square statistic): ");
        double cdf = ChiSquareDistribution::cdf(x, k);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double a = get_double_input("Enter lower bound: ");
        double b = get_double_input("Enter upper bound: ");
        double area = ChiSquareDistribution::area_between(a, b, k);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << a << " <= X <= " << b << ") = " << area << "\n";
    }
}

// Function to handle Exponential distribution
void handle_exponential() {
    display_header(ExponentialDistribution::get_name());
    std::cout << ExponentialDistribution::get_description() << "\n\n";

    double lambda = get_double_input("Enter λ (lambda - rate parameter, e.g., 0.5 events per unit time): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x (time): ");
        double pdf = ExponentialDistribution::pdf(x, lambda);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x (time): ");
        double cdf = ExponentialDistribution::cdf(x, lambda);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double a = get_double_input("Enter lower bound: ");
        double b = get_double_input("Enter upper bound: ");
        double area = ExponentialDistribution::area_between(a, b, lambda);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << a << " <= X <= " << b << ") = " << area << "\n";
    }
}

// Function to handle Beta distribution
void handle_beta() {
    display_header(BetaDistribution::get_name());
    std::cout << BetaDistribution::get_description() << "\n\n";

    double alpha = get_double_input("Enter α (alpha > 0): ");
    double beta_param = get_double_input("Enter β (beta > 0): ");

    std::cout << "\nChoose calculation type:\n";
    std::cout << "  1. PDF at a specific point\n";
    std::cout << "  2. CDF (probability X <= x)\n";
    std::cout << "  3. Area between two points\n";
    int choice = get_int_input("Enter choice: ");

    std::cout << std::fixed << std::setprecision(8);

    if (choice == 1) {
        double x = get_double_input("Enter x (between 0 and 1): ");
        double pdf = BetaDistribution::pdf(x, alpha, beta_param);
        std::cout << "\nResult:\n";
        std::cout << "  PDF at x = " << x << " is " << pdf << "\n";
    } else if (choice == 2) {
        double x = get_double_input("Enter x (between 0 and 1): ");
        double cdf = BetaDistribution::cdf(x, alpha, beta_param);
        std::cout << "\nResult:\n";
        std::cout << "  P(X <= " << x << ") = " << cdf << "\n";
    } else if (choice == 3) {
        double a = get_double_input("Enter lower bound (between 0 and 1): ");
        double b = get_double_input("Enter upper bound (between 0 and 1): ");
        double area = BetaDistribution::area_between(a, b, alpha, beta_param);
        std::cout << "\nResult:\n";
        std::cout << "  P(" << a << " <= X <= " << b << ") = " << area << "\n";
    }
}

// Function to handle combination
void handle_combination() {
    display_header(CombinationCalculator::get_name());
    std::cout << CombinationCalculator::get_description() << "\n\n";

    int n = get_int_input("Enter n (total number of items): ");
    int r = get_int_input("Enter r (number of items to choose): ");

    double result = CombinationCalculator::calculate(n, r);

    std::cout << std::fixed << std::setprecision(0);
    std::cout << "\nResult:\n";
    std::cout << "  C(" << n << ", " << r << ") = " << result << "\n";
    std::cout << "  (Number of ways to choose " << r << " items from " << n << " items)\n";
}

// Function to handle permutation
void handle_permutation() {
    display_header(PermutationCalculator::get_name());
    std::cout << PermutationCalculator::get_description() << "\n\n";

    int n = get_int_input("Enter n (total number of items): ");
    int r = get_int_input("Enter r (number of items to arrange): ");

    double result = PermutationCalculator::calculate(n, r);

    std::cout << std::fixed << std::setprecision(0);
    std::cout << "\nResult:\n";
    std::cout << "  P(" << n << ", " << r << ") = " << result << "\n";
    std::cout << "  (Number of ways to arrange " << r << " items from " << n << " items)\n";
}

// Function to handle factorial
void handle_factorial() {
    display_header(FactorialCalculator::get_name());
    std::cout << FactorialCalculator::get_description() << "\n\n";

    int n = get_int_input("Enter n: ");

    double result = FactorialCalculator::calculate(n);

    std::cout << std::fixed << std::setprecision(0);
    std::cout << "\nResult:\n";
    std::cout << "  " << n << "! = " << result << "\n";
}

// Main menu
void display_main_menu() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════╗\n";
    std::cout << "║   PROBABILITY DISTRIBUTION CALCULATOR      ║\n";
    std::cout << "╚════════════════════════════════════════════╝\n";
    std::cout << "\n";
    std::cout << "═══ DISCRETE DISTRIBUTIONS ═══\n";
    std::cout << "  1. Binomial Distribution\n";
    std::cout << "  2. Negative Binomial Distribution\n";
    std::cout << "  3. Poisson Distribution\n";
    std::cout << "  4. Hypergeometric Distribution\n";
    std::cout << "  5. Geometric Distribution\n";
    std::cout << "\n";
    std::cout << "═══ CONTINUOUS DISTRIBUTIONS ═══\n";
    std::cout << "  6. Normal (Gaussian) Distribution\n";
    std::cout << "  7. Uniform Distribution\n";
    std::cout << "  8. Student's t-Distribution\n";
    std::cout << "  9. F-Distribution\n";
    std::cout << " 10. Gamma Distribution\n";
    std::cout << " 11. Chi-Square Distribution\n";
    std::cout << " 12. Exponential Distribution\n";
    std::cout << " 13. Beta Distribution\n";
    std::cout << "\n";
    std::cout << "═══ SPECIAL FUNCTIONS ═══\n";
    std::cout << " 14. Combination (nCr)\n";
    std::cout << " 15. Permutation (nPr)\n";
    std::cout << " 16. Factorial (n!)\n";
    std::cout << "\n";
    std::cout << " 0. Exit\n";
    std::cout << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    while (true) {
        display_main_menu();
        int choice = get_int_input("Enter your choice: ");

        switch (choice) {
            case 1:
                handle_binomial();
                break;
            case 2:
                handle_negative_binomial();
                break;
            case 3:
                handle_poisson();
                break;
            case 4:
                handle_hypergeometric();
                break;
            case 5:
                handle_geometric();
                break;
            case 6:
                handle_normal();
                break;
            case 7:
                handle_uniform();
                break;
            case 8:
                handle_student_t();
                break;
            case 9:
                handle_f_distribution();
                break;
            case 10:
                handle_gamma();
                break;
            case 11:
                handle_chi_square();
                break;
            case 12:
                handle_exponential();
                break;
            case 13:
                handle_beta();
                break;
            case 14:
                handle_combination();
                break;
            case 15:
                handle_permutation();
                break;
            case 16:
                handle_factorial();
                break;
            case 0:
                std::cout << "\nThank you for using the Probability Calculator!\n";
                return 0;
            default:
                std::cout << "\nInvalid choice. Please try again.\n";
        }

        std::cout << "\nPress Enter to continue...";
        std::cin.get();
    }

    return 0;
}
