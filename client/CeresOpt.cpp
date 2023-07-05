#include "CeresOpt.h"

bool CeresCostFunction::Evaluate(const double *parameters, double *cost,
                                 double *gradient) const {
  const Eigen::VectorXd x =
      Eigen::Map<const Eigen::VectorXd>(parameters, m_objf->degreesOfFreedom());
  m_objf->setPositions(x);
  *cost = m_objf->getEnergy();

  if (gradient != nullptr) {
    Eigen::VectorXd grad = m_objf->getGradient();
    std::copy(grad.data(), grad.data() + grad.size(), gradient);
  }

  return true;
}
int CeresSolver::step(double a_maxMove) {
  m_options.max_num_iterations = 1; // for a single step
  ceres::GradientProblemSolver::Summary summary;
  ceres::Manifold *manifold_ptr = nullptr;
  ceres::GradientProblem problem(m_problem.get(), manifold_ptr);
  ceres::Solve(m_options, problem, nullptr, &summary);
  return summary.termination_type == ceres::CONVERGENCE ? 0 : 1;
}

int CeresSolver::run(size_t a_maxIterations, double a_maxMove) {
  m_options.max_num_iterations = a_maxIterations;
  ceres::GradientProblemSolver::Summary summary;
  ceres::Manifold *manifold_ptr = nullptr;
  ceres::GradientProblem problem(m_problem.get(), manifold_ptr);
  ceres::Solve(m_options, problem, nullptr, &summary);
  return summary.termination_type == ceres::CONVERGENCE ? 0 : 1;
}
