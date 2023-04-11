#include "iwo.h"
#include <random>
#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>

using ms = std::chrono::milliseconds;

void printPopulation(Population& weeds) {
  for (auto& weed : weeds) {
    for (auto& gene : weed)
      std::cout << ' ' << gene;

    std::cout << "; f(x) = " << weed.fitness() << '\n';
  }
}

double sphere(double* values, unsigned size) {
  double sum = 0;
  for (int i = 0; i < size; ++i)
    sum += values[i] * values[i];

  return sum;
}

std::random_device rd;
std::mt19937 gen(rd());

double testMaxSeeds(unsigned dimension,
                    unsigned maxSeeds,
                    double min,
                    double max,
                    unsigned& iterations) {
  std::vector<Interval> limits;
  for (int i = 0; i < dimension; ++i)
    limits.emplace_back(min, max);

  Population weeds(dimension, sphere, 12, limits);

  const unsigned MaxIterations = 16384;
  const unsigned StagnationThreshold = 1024;
  const unsigned NLModIndex = 2;
  const double Eps = 1e-6;
  const double maxSpread = 6.;
  const double minSpread = 1e-10;

  bool stagnationHappened = false;
  unsigned stagnationCount = 0;
  double best = 0, prevBest = std::numeric_limits<double>::max();
  unsigned i;

  for (i = 0; i < MaxIterations && !stagnationHappened; ++i) {
    double spread =
        pow((double)(MaxIterations - i) / MaxIterations, NLModIndex) *
            (maxSpread - minSpread) +
        minSpread;

    weeds.seed(1, 6);
    weeds.spread([spread]() {
      return std::normal_distribution<double>(0, spread)(gen);
    });
    weeds.keep_fittest(maxSeeds);

    prevBest = best;
    best = weeds.best().second;

    if (fabs(best) < Eps)
      stagnationHappened = true;
  }

  auto solution = weeds.best();

  iterations = i;
  /*
  std::cout << "Found solution in " << i << " iterations.\n";
  std::cout << "The found solution is:";
  for (double& gene : solution.first)
    std::cout << ' ' << gene;
  std::cout << "\nThe final fitness function value is: " << solution.second
            << '\n';
  */

  return 0;
}

void runTests(unsigned numTests, unsigned seedCount, double* sumIterations) {
  IWOOptimizer optimizer;

  for (int i = 0; i < numTests; ++i) {
    unsigned iterations;
    optimizer.setMaxSeeds(seedCount);
    auto result = optimizer.test_optimize(
        sphere, {Interval(-10, 10), Interval(-10, 10)}, 0., 1e-6);
    *sumIterations += std::get<2>(result);
  }
}

int main() {
  const unsigned Dimension = 2;
  const unsigned Tests = 16;
  const unsigned NumThreads = 4;
  const std::vector<unsigned> seedCounts = {2,  4,  8,  16,  24,
                                            32, 48, 64, 128, 256};
  IWOOptimizer MahBoi;
  std::chrono::steady_clock clock;
  std::fstream file("data/numSeeds-avgIters.csv");

  for (unsigned seedCount : seedCounts) {
    double sumIterations = 0;
    std::vector<std::thread> threads(4);
    auto tick = clock.now();

    for (int i = 0; i < NumThreads; ++i)
      threads[i] = std::thread(runTests, Tests, seedCount, &sumIterations);

    for (int i = 0; i < NumThreads; ++i)
      threads[i].join();

    auto tock = clock.now();
    double time = std::chrono::duration_cast<ms>(tock - tick).count() / 1000.;
    std::cout << "Convergence found in approximately "
              << sumIterations / (NumThreads * Tests)
              << " iterations.\n With time required: " << time << " s.\n";
    file << seedCount << ", " << time / (sumIterations / (NumThreads * Tests))
         << '\n';
  }
}