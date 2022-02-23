#ifndef PARSEHELPERS_H
#define PARSEHELPERS_H

#include <optional>

#include "Matter.h"
#include "Parameters.h"
#include "Log.h"

namespace helper_functions {
    /**
     * \brief Parse a string into values
     *
     * @param std::string A thing to be parsed
     */
     template <typename T>
     std::vector<T> get_val_from_string(const std::string &line, std::optional<size_t> nelements = std::nullopt);
    /**
     * \brief Split a string into constituent strings
     *
     * Based on https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
     *
     * @param std::string A thing to be parsed
     */
     std::vector<std::string> get_split_strings(const std::string &line);
    } // namespace helper_functions
#endif /* PARSEHELPERS_H */
