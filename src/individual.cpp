#include "iwo.h"
#include <cstring>
#include <random>
#include <limits>

Individual::Individual(unsigned dim,
                       double (*fitness)(double*, unsigned),
                       const DistributionVector& distributions)
    : m_size(dim),
      m_fitnessfunc(fitness),
      m_values(new double[dim]),
      m_has_changed(true) {
  for (int i = 0; i < dim; ++i)
    m_values[i] = distributions[i]();
}

Individual::Individual(const Individual& other)
    : m_size(other.m_size),
      m_fitness(other.m_fitness),
      m_has_changed(other.m_has_changed),
      m_fitnessfunc(other.m_fitnessfunc),
      m_values(new double[other.m_size]) {
  memcpy(m_values, other.m_values, other.m_size * sizeof(double));
}

Individual::Individual(Individual&& other)
    : m_size(other.m_size),
      m_fitness(other.m_fitness),
      m_fitnessfunc(other.m_fitnessfunc),
      m_has_changed(other.m_has_changed),
      m_values(other.m_values) {
  other.m_values = nullptr;
}

Individual& Individual::operator=(const Individual& other) {
  m_size = other.m_size;
  m_fitness = other.m_fitness;
  m_fitnessfunc = other.m_fitnessfunc;
  m_has_changed = other.m_has_changed;

  delete[] m_values;
  m_values = new double[m_size];
  memcpy(m_values, other.m_values, m_size * sizeof(double));

  return *this;
}

Individual& Individual::operator=(Individual&& other) {
  delete[] m_values;
  m_values = other.m_values;
  other.m_values = nullptr;

  m_size = other.m_size;
  m_fitness = other.m_fitness;
  m_fitnessfunc = other.m_fitnessfunc;
  m_has_changed = other.m_has_changed;

  return *this;
}

Individual::~Individual() {
  delete[] m_values;
}

unsigned Individual::dim() const {
  return m_size;
}

double& Individual::operator[](unsigned index) {
  return m_values[index];
}

double Individual::operator[](unsigned index) const {
  return m_values[index];
}

double Individual::fitness() {
  if (m_has_changed)
    m_fitness = m_fitnessfunc(m_values, m_size);

  m_has_changed = false;
  return m_fitness;
}

double* Individual::begin() const {
  return m_values;
}

double* Individual::end() const {
  return m_values + m_size;
}

void Individual::move(DistributionFunc distribution) {
  for (double& value : *this)
    value += distribution();

  m_has_changed = true;
}
