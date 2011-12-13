//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PREFACTOR_H
#define PREFACTOR_H

#include "Eigen.h"
#include "Matter.h"

#include "Prefactor.h"
#include "Parameters.h"

int getPrefactors(Parameters *parameters, Matter *min1, Matter *saddle, Matter *min2, double &pref1, double &pref2);
VectorXi movedAtoms(Parameters* parameters, Matter *min1, Matter *saddle, Matter *min2);
VectorXd removeZeroFreqs(Parameters *parameters, VectorXd freqs);

#endif
