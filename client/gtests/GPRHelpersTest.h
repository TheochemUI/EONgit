#ifndef GPRHELPERSTEST_H
#define GPRHELPERSTEST_H

#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <functional>
#include <memory>
#include <numeric>

#include "../GPRHelpers.h"
#include "../HelperFunctions.h"
#include "../Job.h"
#include "../Log.h"
#include "../Matter.h"
#include "../MinModeSaddleSearch.h"
#include "../Parameters.h"
#include "../Potential.h"
#include "../external/icecream.hpp"
#include "../potentials/LJ/LJ.h"
#include "../potentials/Morse/Morse.h"

namespace tests {
class GPRHelpersTest : public ::testing::Test {
 public:
  GPRHelpersTest();
  virtual ~GPRHelpersTest();
  std::string reactantFilename;
  std::string productFilename;
  const long double threshold{1e3 * std::numeric_limits<double>::epsilon()};

  std::unique_ptr<Parameters> parameters;
  std::unique_ptr<Matter> initmatter;
  std::unique_ptr<Matter> finalmatter;
  std::function<bool(const gpr::EigenMatrix &lhs, const gpr::EigenMatrix &rhs)>
      comparer;
};
} /* namespace tests */

#endif /* GPRHELPERSTEST_H */
