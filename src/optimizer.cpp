#include "iwo.h"
#include "optimizer.h"
#include <cmath>

IWOOptimizer::IWOOptimizer() : m_random(std::random_device()()) {}

std::tuple<Individual&, double, unsigned> IWOOptimizer::test_optimize(
    double (*objective)(double*, unsigned),
    const std::vector<Interval>& limits,
    double reference,
    double precision) {
  if (m_herbs)
    delete m_herbs;
  m_herbs =
      new Population((unsigned)limits.size(), objective, m_init_size, limits);

  unsigned i;
  bool precisionReached = false;

  for (i = 0; i < m_max_iterations && !precisionReached; ++i) {
    double variance =
        pow((double)(m_max_iterations - i) / m_max_iterations, m_nlmindex) *
            (m_max_variance - m_min_variance) +
        m_min_variance;
    m_herbs->seed(m_min_seeds, m_max_seeds);
    m_herbs->spread([this, variance]() {
      return std::normal_distribution<double>(0., variance)(m_random);
    });
    m_herbs->keep_fittest(m_max_seeds);

    if (fabs(m_herbs->best().second - reference) < precision)
      precisionReached = true;
  }

  auto solution = m_herbs->best();
  return {solution.first, solution.second, i};
}

void IWOOptimizer::setInitPopulationSize(unsigned count) {
  m_init_size = count;
}

void IWOOptimizer::setMinSeeds(unsigned count) {
  m_min_seeds = count;
}

void IWOOptimizer::setMaxSeeds(unsigned count) {
  m_max_seeds = count;
}

void IWOOptimizer::setMinVariance(double value) {
  m_min_variance = value;
}

void IWOOptimizer::setMaxVariance(double value) {
  m_max_variance = value;
}

void IWOOptimizer::setPopulationCap(unsigned count) {
  m_max_size = count;
}

void IWOOptimizer::setMaxIterations(unsigned count) {
  m_max_iterations = count;
}

LFIWOOptimizer::LFIWOOptimizer() : IWOOptimizer() {}
