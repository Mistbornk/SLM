#pragma once
#include <string>
#include <iostream>

inline bool verbose = false;

inline void verbose_log(const std::string& message) {
    if (verbose) {
        std::cout << message << std::endl;
    }
}