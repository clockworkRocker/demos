#ifndef IWO_H
#define IWO_H

#include <functional>
#include <tuple>
#include <vector>
#include <random>

using DistributionFunc = std::function<double()>;
using DistributionVector = std::vector<DistributionFunc>;

#include "individual.h"
#include "population.h"
#include "optimizer.h"

#endif  // IWO_H