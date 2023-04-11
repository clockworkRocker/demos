#ifndef IWO_OPTIMIZER_H
#define IWO_OPTIMIZER_H

class IWOOptimizer {
 public:
  static constexpr double DefaultMinSpread = 1e-6;
  static constexpr double DefaultMaxSpread = 4.;
  enum {
    DefaultInitialPopulationSize = 12,
    DefaultMinSeeds = 1,
    DefaultMaxSeeds = 8,
    DefaultNLMIndex = 2,
    DefaultPopulationCap = 24,
    DefaultMaxIterations = 16384,
  };

 public:
  IWOOptimizer();
  /// @brief Test the optimizer
  /// @param objective An objective function with the known minimum
  /// @param limits Intervals that define the search space
  /// @param reference Reference value of the known objective function minimum
  /// @param precision The precision with which the solution needs to be found
  /// @return The set of the best individual of the popilation, its fitness
  ///         value and iterations used
  std::tuple<Individual&, double, unsigned> test_optimize(
      double (*objective)(double*, unsigned),
      const std::vector<Interval>& limits,
      double reference,
      double precision);

 public:
  void setInitPopulationSize(unsigned count);
  void setMinSeeds(unsigned count);
  void setMaxSeeds(unsigned count);
  void setMinVariance(double value);
  void setMaxVariance(double value);
  void setPopulationCap(unsigned count);
  void setMaxIterations(unsigned count);

 protected:
  unsigned m_init_size = DefaultInitialPopulationSize;
  unsigned m_max_size = DefaultPopulationCap;
  unsigned m_min_seeds = DefaultMinSeeds;
  unsigned m_max_seeds = DefaultMaxSeeds;
  double m_min_variance = DefaultMinSpread;
  double m_max_variance = DefaultMaxSpread;
  double m_nlmindex = DefaultNLMIndex;
  double m_max_iterations = DefaultMaxIterations;
  std::mt19937 m_random;

  Population* m_herbs = nullptr;
};

/* ========================================================================= */

class LFIWOOptimizer : public IWOOptimizer {
 public:
  LFIWOOptimizer();

 public:
  /// @brief Test the optimizer
  /// @param objective An objective function with the known minimum
  /// @param limits Intervals that define the search space
  /// @param reference Reference value of the known objective function minimum
  /// @param precision The precision with which the solution needs to be found
  /// @return The set of the best individual of the popilation, its fitness
  ///         value and iterations used
  std::tuple<Individual&, double, unsigned> test_optimize(
      double (*objective)(double*, unsigned),
      const std::vector<Interval>& limits,
      double reference,
      double precision);
};

#endif  // IWO_OPTIMIZER_H