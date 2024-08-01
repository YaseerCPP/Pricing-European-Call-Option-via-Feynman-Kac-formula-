#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

// Function to calculate the payoff of a European call option
double callPayoff(double S, double K) {
    return std::max(S - K, 0.0);
}

// Function to simulate the stock price at maturity
double simulateStockPrice(double S0, double r, double sigma, double T, std::mt19937 &gen, std::normal_distribution<> &d) {
    double ST = S0 * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * d(gen));
    return ST;
}

// Monte Carlo simulation to price a European call option
double monteCarloCallPrice(double S0, double K, double r, double sigma, double T, int numSimulations) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);

    double sumPayoff = 0.0;

    for (int i = 0; i < numSimulations; ++i) {
        double ST = simulateStockPrice(S0, r, sigma, T, gen, d);
        double payoff = callPayoff(ST, K);
        sumPayoff += payoff;
    }

    double discountedPayoff = (sumPayoff / numSimulations) * std::exp(-r * T);
    return discountedPayoff;
}

int main() {
    // Option parameters
    double S0 = 100.0;   // Initial stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free interest rate
    double sigma = 0.2;  // Volatility
    double T = 1.0;      // Time to maturity in years
    int numSimulations = 100000; // Number of Monte Carlo simulations

    // Calculate the option price
    double callPrice = monteCarloCallPrice(S0, K, r, sigma, T, numSimulations);

    std::cout << "European Call Option Price: " << callPrice << std::endl;

    return 0;
}
