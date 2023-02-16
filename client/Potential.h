#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Eigen.h"
#include "Parameters.h"
#include <algorithm>
#include <limits>
#include <memory>

class Potential {
private:
  // should be const
  PotType ptype;
  Parameters *params;

public:
  Potential(Parameters *parameters) : ptype{parameters->potential} {}
  ~Potential() = default;

  static int fcalls;
  static int fcallsTotal;
  static int wu_fcallsTotal;
  static double totalUserTime;

// Does not take into account the fixed / free atoms
  std::pair<double, AtomMatrix> virtual get_ef(const AtomMatrix pos,
                                               const VectorXi atmnrs,
                                               const Matrix3d box) = 0;
  PotType getType() { return this->ptype; };
};

namespace helper_functions {
Potential *makePotential(Parameters *params);
} // namespace helper_functions

// int Potential::fcalls = 0;
// int Potential::fcallsTotal = 0;
// int Potential::wu_fcallsTotal = 0;
// double Potential::totalUserTime = 0.0;
#endif
