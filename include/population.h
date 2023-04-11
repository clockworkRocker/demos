#ifndef IWO_POPULATION_H
#define IWO_POPULATION_H

struct byFitness {
  bool operator()(Individual& lhs, Individual& rhs) const {
    return lhs.fitness() < rhs.fitness();
  }
};

struct Interval {
  double a;
  double b;

  explicit Interval(double _a = 0, double _b = 0);
};

class Population {
 public:  // * Constructors and stuff
  Population(unsigned dim,
             double (*fitnessfunc)(double*, unsigned),
             unsigned size,
             const std::vector<Interval>& limits);

 public:  // * Access methods
  std::pair<Individual&, double> best();
  double average_fitness();

 public:  // * Iterators
  typedef std::vector<Individual>::iterator herb_iterator;
  herb_iterator begin();
  herb_iterator end();

 public:  // * Optimization steps
  /// @brief Create new seeds
  void seed(unsigned min_seeds, unsigned max_seeds);

  /// @brief Use the distribution function to spread the seeds of the population
  void spread(DistributionFunc distribution);

  /// @brief Limit the population to the new maximal size keeping the fittest
  ///        individuals
  void keep_fittest(unsigned count);

 private:
  bool in_limits(const Individual& individual);

 private:
  double (*m_fitnessfunc)(double*, unsigned);
  unsigned m_num_herbs;

  std::vector<Interval> m_limits;
  std::vector<DistributionFunc> m_searchspace;
  std::vector<Individual> m_herbs;
};

#endif  // IWO_POPULATION_H