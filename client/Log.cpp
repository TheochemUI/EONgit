#include "Log.h"

void log(const std::string& message, spdlog::level::level_enum level) {
    auto consoleLogger = spdlog::get("console");
    if (!consoleLogger) {
        consoleLogger = spdlog::stdout_color_mt("console");
    }
    consoleLogger->log(level, message);
}

void log_file(const std::string& message, spdlog::level::level_enum level) {
    auto fileLogger = spdlog::get("file");
    if (!fileLogger) {
        fileLogger = spdlog::basic_logger_mt("file", "log.txt");
    }
    fileLogger->log(level, message);
}
