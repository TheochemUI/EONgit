/*
 * GPRHelpersTest.cpp
 *
 *  Created on: 07 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../potentials/LJ/LJ.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Job.h"
#include "../Log.h"
#include "../Matter.h"
#include "../MinModeSaddleSearch.h"
#include "../Parameters.h"
#include "GPRHelpersTest.h"

namespace tests {

GPRHelpersTest::GPRHelpersTest() {
  // TODO Auto-generated constructor stub
}

GPRHelpersTest::~GPRHelpersTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRHelpersTest, TestMatter) {
  string confile("pos.con");
  double energy_morse_ref{-1775.79}; // matter->getPotentialEnergy()
  Parameters *pmorse = new Parameters;
  pmorse->potential = "morse_pt";
  pmorse->LogPotential = false;
  Matter *pm = new Matter(pmorse);
  pm->con2matter(confile);
  Morse pot_morse;
  pot_morse.setParams(pmorse); // *mostly* pointless
  EXPECT_EQ(pot_morse.getName(), "morse_pt"s)
    << "Potential name is incorrect [Morse]";
  double energy_morse{0};
  AtomMatrix forces_morse = AtomMatrix::Constant(pm->numberOfAtoms(), 3, 0);
  auto egf_morse = helper_functions::energy_and_forces(pm, &pot_morse);
  energy_morse = std::get<double>(egf_morse);
  forces_morse = std::get<AtomMatrix>(egf_morse);
  EXPECT_EQ(forces_morse, pm->getForces())
      << "Forces do not match [Morse]";
  EXPECT_EQ(energy_morse, pm->getPotentialEnergy())
      << "Potential energy does not match [Morse]";
  Parameters *plj = new Parameters;
  plj->potential = "lj";
  plj->LogPotential = false;
  Matter *pljm = new Matter(plj);
  pljm->con2matter(confile);
  LJ ljpot;
  ljpot.setParams(plj);
  EXPECT_EQ(ljpot.getName(), "lj"s)
    << "Potential name is incorrect [LJ]";
  double energy_lj{0};
  AtomMatrix forces_lj = AtomMatrix::Constant(pljm->numberOfAtoms(), 3, 0);
  auto egf_lj = helper_functions::energy_and_forces(pm, &ljpot);
  energy_lj = std::get<double>(egf_lj);
  forces_lj = std::get<AtomMatrix>(egf_lj);
  pljm->setPotential(&ljpot);
  EXPECT_EQ(pljm->getPotential()->getName(), "lj"s)
    << "Potential name is incorrect in matter [LJ]";
  EXPECT_EQ(forces_lj, pljm->getForces())
      << "Forces do not match [LJ]";
  EXPECT_EQ(energy_lj, pljm->getPotentialEnergy())
      << "Potential energy does not match [LJ]";
  delete pljm;
  delete plj;
  delete pm;
  delete pmorse;
}

} /* namespace tests */
