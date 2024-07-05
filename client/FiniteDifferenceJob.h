#pragma once
#include "Job.h"
#include "Parameters.h"

class FiniteDifferenceJob : public Job {
public:
  FiniteDifferenceJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~FiniteDifferenceJob(void) = default;
  std::vector<std::string> run(void);
};
