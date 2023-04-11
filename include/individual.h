#ifndef IWO_INDIVIDUAL_H
#define IWO_INDIVIDUAL_H

class Individual {
 public:  // Constructors and assignment
  /// @brief Create a new individual and fill it with random numbers
  /// @param dim Input dimension
  /// @param fitness Pointer to the fitness function
  /// @param distributions The distributions to generate the numbers from
  Individual(unsigned dim = 0,
             double (*fitness)(double*, unsigned) = nullptr,
             const DistributionVector& distributions = {});

  /// @brief Copy an individual
  Individual(const Individual& other);

  /// @brief Move an individual object
  Individual(Individual&& other);

  /// @brief Assignment
  Individual& operator=(const Individual& other);

  /// @brief Move-assignment
  Individual& operator=(Individual&& other);

  ~Individual();

 public:  // Access methods
  /// @return Input dimension
  unsigned dim() const;

  /// @return A reference to the element of individual vector with the given
  ///         index
  double& operator[](unsigned index);

  /// @return The value at a given index in the individual vector
  double operator[](unsigned index) const;

  /// @return The value of the fitness function
  double fitness();

 public:  // Iterators
  double* begin() const;
  double* end() const;

 public:  // IWO Methods
  /// @brief Move the individual in the search space by a random distance
  ///         defined by normal distribution with null mean and variance
  void move(DistributionFunc distribution);

 protected:
  unsigned m_size;
  double* m_values;
  double m_fitness;
  double (*m_fitnessfunc)(double*, unsigned);

  bool m_has_changed = false;
};

#endif  // IWO_INDIVIDUAL_H