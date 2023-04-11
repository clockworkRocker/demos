#include "iwo.h"
#include <random>
#include <cassert>
#include <algorithm>

Interval::Interval(double _a, double _b) : a(_a), b(_b) {}

/* ========================================================================= */

Population::Population(unsigned dim,
                       double (*fitnessfunc)(double*, unsigned),
                       unsigned size,
                       const std::vector<Interval>& limits)
    : m_num_herbs(size),
      m_fitnessfunc(fitnessfunc),
      m_herbs(),
      m_searchspace(),
      m_limits(limits) {
  assert(dim == limits.size());

  std::random_device rd;
  std::mt19937 gen(rd());

  for (const Interval& limit : limits)
    m_searchspace.push_back([limit, &gen]() {
      return std::uniform_real_distribution<double>(limit.a, limit.b)(gen);
    });

  for (int i = 0; i < size; ++i)
    m_herbs.emplace_back(dim, fitnessfunc, m_searchspace);
}

std::pair<Individual&, double> Population::best() {
  double best = m_herbs[0].fitness();
  Individual& bestHerb = m_herbs[0];

  for (Individual& herb : m_herbs)
    if (herb.fitness() < best) {
      best = herb.fitness();
      bestHerb = herb;
    }

  return {bestHerb, best};
}

double Population::average_fitness() {
  double sum = 0;

  for (auto& herb : m_herbs) {
    sum += herb.fitness();
  }

  return sum / m_herbs.size();
}

Population::herb_iterator Population::begin() {
  return m_herbs.begin();
}

Population::herb_iterator Population::end() {
  return m_herbs.end();
}

void Population::seed(unsigned min_seeds, unsigned max_seeds) {
  std::vector<Individual> withSeeds(m_herbs);
  unsigned num_seeds;
  auto tops = std::minmax_element(m_herbs.begin(), m_herbs.end(), byFitness());
  double bestFit = tops.first->fitness();
  double worstFit = tops.second->fitness();

  for (auto& herb : m_herbs) {
    num_seeds = min_seeds + (max_seeds - min_seeds) / (worstFit - bestFit) *
                                (worstFit - herb.fitness());
    for (int j = 0; j < num_seeds; ++j)
      withSeeds.emplace_back(herb);
  }

  m_herbs = std::move(withSeeds);
}

void Population::spread(DistributionFunc distribution) {
  for (auto iter = m_herbs.begin() + m_num_herbs; iter != m_herbs.end(); ++iter)
    iter->move(distribution);
}

void Population::keep_fittest(unsigned count) {
  m_herbs.erase(
      std::remove_if(m_herbs.begin(), m_herbs.end(),
                     [this](const Individual& i) { return !in_limits(i); }),
      m_herbs.end());

  std::sort(m_herbs.begin(), m_herbs.end(), byFitness());
  if (m_herbs.size() > count)
    m_herbs.resize(count);

  m_num_herbs = std::min(m_herbs.size(), (size_t)count);
}

bool Population::in_limits(const Individual& individual) {
  for (int i = 0; i < individual.dim(); ++i)
    if (individual[i] < m_limits[i].a || individual[i] > m_limits[i].b)
      return false;

  return true;
}
