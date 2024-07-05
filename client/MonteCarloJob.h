#pragma once

#include "Job.h"
#include "Parameters.h"

class MonteCarloJob : public Job {
public:
  MonteCarloJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    log = spdlog::get("combi");
  }
  ~MonteCarloJob(void) = default;
  std::vector<std::string> run(void);

private:
  std::shared_ptr<spdlog::logger> log;
};
