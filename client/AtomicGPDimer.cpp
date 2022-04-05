// an interface to the gpdimer library

#include "AtomicGPDimer.h"
#include "HelperFunctions.h"
#include "GPRHelpers.h"
#include "Log.h"
#include <cassert>
#include <cmath>

#include "subprojects/gprdimer/gpr/AtomicDimer.h"
#include "subprojects/gprdimer/gpr/Enums.h"
#include "subprojects/gprdimer/structures/Structures.h"
#include "subprojects/gprdimer/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gprdimer/gpr/covariance_functions/ConstantCF.h"
#include "subprojects/gprdimer/gpr/covariance_functions/SexpatCF.h"
#include "subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

AtomicGPDimer::AtomicGPDimer(Matter *matter, Parameters *params) {
  parameters = params;
  matterCenter = new Matter(parameters);
  *matterCenter = *matter;
  p = helper_functions::eon_parameters_to_gpr(params);
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = matter->getCell()(i);
  }
}

AtomicGPDimer::~AtomicGPDimer() { delete matterCenter; }

void AtomicGPDimer::compute(Matter *matter,
                            AtomMatrix initialDirectionAtomMatrix) {
  *matterCenter = *matter;
  auto config_data = helper_functions::eon_matter_to_frozen_conf_info(matterCenter,
                                                             parameters->gprActiveRadius);
  auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
  auto R_init = std::get<gpr::Coord>(config_data);
  // R_init.resize(1, matterCenter->getPositionsFree().size());
  // R_init.assignFromEigenMatrix(matterCenter->getPositionsFreeV());
  init_middle_point.clear();
  init_middle_point.R = R_init;
  init_observations.clear();
  orient_init.clear();
  orient_init.resize(matterCenter->getPositionsFree().rows(),
                     matterCenter->getPositionsFree().cols());
    AtomMatrix freeOrient(matterCenter->numberOfFreeAtoms(),3);
    int i,j = 0;
    for(i=0; i<matterCenter->numberOfAtoms(); i++)
    {
        if(!matterCenter->getFixed(i))
        {
           freeOrient.row(j) = initialDirectionAtomMatrix.row(i);
            j++;
            if (j==matterCenter->numberOfFixedAtoms()){
              break;
            }
        }
    }
  // orient_init.assignFromEigenMatrix(freeOrient);
orient_init.resize(1, freeOrient.rows() *
                freeOrient.cols());
  for(int i,counter = 0; i < freeOrient.rows(); ++i) {
    for(int j = 0; j < freeOrient.cols(); ++j) {
      orient_init[counter++] = freeOrient(i, j);
    }
  }
  atomic_dimer.initialize(p, init_observations, init_middle_point, orient_init,
                          atoms_config);

  Potential *potential= Potential::getPotential(parameters);
  atomic_dimer.execute(*potential);
  // Forcefully set the right positions
  matter->setPositionsFreeV(atomic_dimer.getFinalCoordOfMidPoint());
  this->totalIterations=atomic_dimer.getIterations();
  this->totalForceCalls=atomic_dimer.getTotalForceCalls();
  return;
}

double AtomicGPDimer::getEigenvalue() {
  return atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector() {
  return atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
