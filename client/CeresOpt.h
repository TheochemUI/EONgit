#ifndef CERESOPT_H_
#define CERESOPT_H_

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"
#include <ceres/ceres.h>

class CeresCostFunction : public ceres::FirstOrderFunction {
public:
  CeresCostFunction(std::shared_ptr<ObjectiveFunction> a_objf)
      : m_objf{a_objf} {}
  virtual ~CeresCostFunction() {}
  bool Evaluate(const double *parameters, double *cost,
                double *gradient) const override;
  int NumParameters() const override { return m_objf->degreesOfFreedom(); }

private:
  std::shared_ptr<ObjectiveFunction> m_objf;
};

class CeresSolver final : public Optimizer {
public:
  CeresSolver(std::shared_ptr<ObjectiveFunction> a_objf,
              std::shared_ptr<Parameters> a_params)
      : Optimizer(a_objf, OptType::ceres, a_params),
        m_problem(std::make_unique<CeresCostFunction>(a_objf)) {
    if (spdlog::get("ceres")) {
      m_log = spdlog::get("ceres");
    } else {
      m_log = spdlog::basic_logger_st("ceres", "_ceres.log", true);
    }
    m_log->set_pattern("[%l] [CERES] %v");

    m_options.minimizer_progress_to_stdout = true;
    m_options.max_num_iterations =
        a_params->optMaxIterations; // Set your options here
    // m_options.linear_solver_type = ceres::DENSE_QR;
  }

  ~CeresSolver() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  std::unique_ptr<CeresCostFunction> m_problem;
  ceres::GradientProblemSolver::Options m_options;
  std::shared_ptr<spdlog::logger> m_log;
};

#endif // CERESOPT_H_
