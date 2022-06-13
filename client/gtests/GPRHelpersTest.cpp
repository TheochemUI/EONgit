/*
 * GPRHelpersTest.cpp
 *
 *  Created on: 07 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "GPRHelpersTest.h"

namespace tests {

GPRHelpersTest::GPRHelpersTest() {
  reactantFilename = helper_functions::getRelevantFile("reactant.con");
  productFilename = helper_functions::getRelevantFile("product.con");

  parameters = std::make_unique<Parameters>();
  parameters->potential = "morse_pt";
  parameters->nebImages = 7;
  parameters->LogPotential = false;
  log_init(parameters.get(), (char *)"test.log");

  initmatter = std::make_unique<Matter>(parameters.get());
  finalmatter = std::make_unique<Matter>(parameters.get());

  initmatter->con2matter(reactantFilename);
  finalmatter->con2matter(productFilename);
  // Comparer
  this->comparer = [&](const gpr::EigenMatrix &lhs,
                       const gpr::EigenMatrix &rhs) -> bool {
    return lhs.isApprox(rhs, this->threshold);
  };
};

GPRHelpersTest::~GPRHelpersTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRHelpersTest, TestNEBInitPath) {
  // Setup the run
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  ASSERT_PRED2(comparer, imgArray.front().getForces(),
               this->initmatter->getForces())
      << "Energies of the first image don't match";
  ASSERT_PRED2(comparer, imgArray.back().getForces(),
               this->finalmatter->getForces())
      << "Energies of the last image don't match";
}

// TEST_F(GPRHelpersTest, TestNEBInitObs) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto obspath = helper_functions::prepInitialObs(imgArray);
//   // Assertions
//   const int nimages = this->parameters->nebImages;
//   const int totImages = nimages + 2; // Final and end
//   EXPECT_EQ(obspath.E.getSize(), totImages)
//       << "Number of elements of energies do not match";
//   EXPECT_EQ(obspath.R.getNumPoints(), totImages)
//       << "Number of elements of R matrices do not match";
//   EXPECT_EQ(obspath.G.getNumPoints(), totImages)
//       << "Number of elements of G matrices do not match";
//   // obspath.printSizes();
// }

// TEST_F(GPRHelpersTest, TestPathLength) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   double plen = helper_functions::get_path_length(imgArray);
//   double true_plen = imgArray.front().distanceTo(imgArray.back());
//   EXPECT_NEAR(plen, true_plen, this->threshold)
//       << "Path length is incorrect";
// }

// TEST_F(GPRHelpersTest, TestUnravelFreeV) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto ufreecoords = helper_functions::unravel_free_coordsV(imgArray);
//   const size_t nfree = imgArray.front().numberOfFreeAtoms();
//   const size_t nimgs = imgArray.size();
//   EXPECT_EQ(ufreecoords.size(), nfree*nimgs*3)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   for (size_t idx{0}; idx < nimgs; idx++){
//     EXPECT_EQ(ufreecoords.segment(idx * nfree * 3, nfree * 3),
//     imgArray[idx].getPositionsFreeV())
//       << "Point isn't the same";
//   }
// }

// TEST_F(GPRHelpersTest, TestUnravelV) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto ucoords = helper_functions::unravel_coordsV(imgArray);
//   const size_t npoints = imgArray.front().numberOfAtoms();
//   const size_t nimgs = imgArray.size();
//   EXPECT_EQ(ucoords.size(), npoints*nimgs*3)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   for (size_t idx{0}; idx < nimgs; idx++){
//     EXPECT_EQ(ucoords.segment(idx * npoints * 3, npoints * 3),
//     imgArray[idx].getPositionsV())
//       << "Point isn't the same";
//   }
// }

// TEST_F(GPRHelpersTest, TestCoordGeneration) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   Matter oneThing = imgArray.front();
//   gpr::Coord oneCoord = helper_functions::single_img(oneThing);
//   for (gpr::Index_t n = 0; n < oneCoord.getSize(); ++n) {
//     EXPECT_EQ(oneCoord.getData()[n], oneThing.getPositionsFreeV().data()[n])
//     << "The field has wrong values.";
//   }
// }

// TEST_F(GPRHelpersTest, TestCoordPath) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   gpr::Coord coordPath = helper_functions::prev_path(imgArray);
//   // for (gpr::Index_t n = 0; n < oneCoord.getSize(); ++n) {
//   //   EXPECT_EQ(oneCoord.getData()[n],
//   oneThing.getPositionsFreeV().data()[n]) << "The field has wrong values.";
//   // }
// }

// TEST_F(GPRHelpersTest, TestUnravelFree) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto ufreecoords = helper_functions::unravel_free_coords(imgArray);
//   const size_t nfree = imgArray.front().numberOfFreeAtoms();
//   const size_t nimgs = imgArray.size();
//   size_t nimgrows{0};
//   for (auto img:imgArray){
//     nimgrows += img.numberOfFreeAtoms();
//   }
//   // IC(ufreecoords.rows(), ufreecoords.cols(), ufreecoords.size());
//   EXPECT_EQ(ufreecoords.rows(), nimgrows)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   EXPECT_EQ(ufreecoords.cols(), nfree*3)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   for (size_t idx{0}; idx < nimgs; idx++){
//     EXPECT_EQ(ufreecoords.row(idx), imgArray[idx].getPositionsFree())
//       << "Point isn't the same";
//     // IC(ufreecoords.row(idx), imgArray[idx].getPositionsFree());
//   }
// }

// TEST_F(GPRHelpersTest, TestUnravel) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto ucoords = helper_functions::unravel_coords(imgArray);
//   const size_t npositions = imgArray.front().numberOfAtoms();
//   const size_t nimgs = imgArray.size();
//   size_t nimgrows{0};
//   // for (auto img:imgArray){
//   //   nimgrows += img.numberOfAtoms();
//   // }
//   IC(ucoords.rows(), ucoords.cols(), ucoords.size());
//   // EXPECT_EQ(ucoords.rows(), nimgrows)
//   //   << "Unraveled coordinates have the wrong number of elements!\n";
//   EXPECT_EQ(ucoords.cols(), npositions*3)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   for (size_t idx{0}; idx < nimgs; idx++){
//     EXPECT_EQ(ucoords.row(idx), imgArray[idx].getPositions())
//       << "Point isn't the same";
//     // IC(ucoords.row(idx), imgArray[idx].getPositions());
//   }
// }

} /* namespace tests */
