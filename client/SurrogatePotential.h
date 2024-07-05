#pragma once

#include "Potential.h"

class SurrogatePotential : public Potential {

public:
  SurrogatePotential(PotType a_ptype, std::shared_ptr<Parameters> a_params)
      : Potential(a_ptype, a_params) {}
  virtual ~SurrogatePotential() = default;
  std::tuple<double, AtomMatrix, double> // energy, forces, energy variance
  get_ef_var(const AtomMatrix pos, const Vector<int> atmnrs,
             const Matrix3S box);
  virtual void train_optimize(MatrixType a_features, MatrixType a_targets) = 0;
};
