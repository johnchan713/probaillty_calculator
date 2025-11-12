#pragma once

#include "math_utils.hpp"
#include <string>

namespace ProbCalc {

// Combination (nCr) Calculator
class CombinationCalculator {
public:
    static std::string get_name() { return "Combination (nCr)"; }

    static std::string get_description() {
        return "Calculates the number of ways to choose r items from n items.\n"
               "Order does NOT matter.\n"
               "Formula: nCr = n! / (r! * (n-r)!)\n"
               "Example: Choosing 3 cards from a deck of 52.\n"
               "Example: Selecting 5 people from a group of 10.\n"
               "Example: Lottery combinations.";
    }

    static double calculate(int n, int r) {
        return combination(n, r);
    }
};

// Permutation (nPr) Calculator
class PermutationCalculator {
public:
    static std::string get_name() { return "Permutation (nPr)"; }

    static std::string get_description() {
        return "Calculates the number of ways to arrange r items from n items.\n"
               "Order DOES matter.\n"
               "Formula: nPr = n! / (n-r)!\n"
               "Example: Arranging 3 books from a shelf of 10.\n"
               "Example: Race finishing positions (1st, 2nd, 3rd from n runners).\n"
               "Example: Password permutations.";
    }

    static double calculate(int n, int r) {
        return permutation(n, r);
    }
};

// Factorial Calculator
class FactorialCalculator {
public:
    static std::string get_name() { return "Factorial (n!)"; }

    static std::string get_description() {
        return "Calculates n! = n × (n-1) × (n-2) × ... × 2 × 1\n"
               "Formula: n! = n × (n-1)!\n"
               "Example: Number of ways to arrange n items.\n"
               "Example: 5! = 5 × 4 × 3 × 2 × 1 = 120\n"
               "Note: 0! = 1 by definition.";
    }

    static double calculate(int n) {
        return factorial(n);
    }
};

} // namespace ProbCalc
