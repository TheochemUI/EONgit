#ifndef EON_LOG_H
#define EON_LOG_H

#include "Parameters.h"

// Free function to log to console
void log(const std::string& message, spdlog::level::level_enum level = spdlog::level::debug);
// Free function to log to file
void log_file(const std::string& message, spdlog::level::level_enum level = spdlog::level::debug);

#endif
