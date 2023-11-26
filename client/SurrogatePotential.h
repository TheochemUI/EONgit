#pragma once

#include "Matter.h"
#include "Potential.h"

class SurrogatePotential : public Potential {

public:
  bool failedOptim;
  SurrogatePotential(PotType a_ptype, std::shared_ptr<Parameters> a_params)
    : Potential(a_ptype, a_params), failedOptim{false} {}
  virtual ~SurrogatePotential() = default;
  virtual void train_optimize(Eigen::MatrixXd a_features,
                              Eigen::MatrixXd a_targets) = 0;
  std::tuple<double, AtomMatrix, double> // energy, forces, energy variance
  get_ef_var(const AtomMatrix pos, const VectorXi atmnrs, const Matrix3d box);
};
