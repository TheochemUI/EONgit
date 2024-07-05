#pragma once
#include "Eigen.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"

// Lanczos method to find the lowest curvature mode
class Lanczos : public LowestEigenmode {

public:
  Lanczos(std::shared_ptr<Matter> matter, std::shared_ptr<Parameters> params,
          std::shared_ptr<Potential> pot);
  ~Lanczos() = default;
  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue();
  AtomMatrix getEigenvector();

private:
  AtomMatrix lowestEv;
  double lowestEw;
  std::shared_ptr<spdlog::logger> log;
};
