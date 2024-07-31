#include "EnvHelpers.hpp"
#include <spdlog/spdlog.h>
namespace eonc::def {

std::string get_value_from_env_or_param(const char *env_variable,
                                        const PDef val_def,
                                        const std::string &warning_message,
                                        const bool is_mandatory) {
  const char *env_value = std::getenv(env_variable);
  auto &&[param_value, default_value] = val_def;
  if (env_value != nullptr) {
    return std::string(env_value);
  }

  if (!param_value.empty()) {
    return param_value;
  }

  if (is_mandatory) {
    throw std::runtime_error(
        "Environment variable " + std::string(env_variable) +
        " is not set and no parameter value provided. Please set it in the "
        "configuration or as an environment variable.\n");
  }

  if (!default_value.empty() && !warning_message.empty()) {
    SPDLOG_WARN(warning_message);
  }

  return default_value;
}

} // namespace eonc::def
