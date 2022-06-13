/*
 * ObsTest.cpp
 *
 *  Created on: 05 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "ObsTest.h"
#include "../external/icecream.hpp"

namespace tests {

ObsTest::ObsTest() {
  // TODO Auto-generated constructor stub
}

ObsTest::~ObsTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(ObsTest, TestMatter) {
  string confile("pos.con");
  Parameters *parameters = new Parameters;
  parameters->potential = "morse_pt";
  Matter *matter = new Matter(parameters);
  matter->con2matter(confile);
  gpr::Observation o = helper_functions::eon_matter_to_init_obs(*matter);
  EXPECT_EQ(o.R.extractEigenMatrix().reshaped<Eigen::RowMajor>(matter->numberOfFreeAtoms(), 3), matter->getPositionsFree()) // Figure out the right approach..
      << "Positions do not match";
  EXPECT_EQ((o.G.extractEigenMatrix()*-1).reshaped<Eigen::RowMajor>(matter->numberOfFreeAtoms(), 3), matter->getForcesFree())
      << "Forces do not match";
  EXPECT_EQ(o.E[0], matter->getPotentialEnergy())
      << "Potential energy does not match";
  delete matter;
  delete parameters;
}

} /* namespace tests */
