#include "zx/Measurement.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <optional>
#include <sstream>
#include <tuple>
#include <windows.h>

Measurement::Measurement(bool enabled):
    enabled(enabled){};

void Measurement::addMeasurement(std::string name, std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end) {
    if (enabled) {
        double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
#pragma omp critical
        {
            times[name] += duration;
        }
    }
}

void Measurement::printMeasurements(std::string group, std::string circuit, int parallel_iterations, int iterations, std::string filename) {
    for (auto t: times) {
        std::cout << t.first << " = " << t.second / 1000000.0 << "[ms]" << std::endl;
    }

    std::ofstream     file(filename, std::ios::app);
    auto              now   = std::chrono::system_clock::now();
    std::time_t       now_c = std::chrono::system_clock::to_time_t(now);
    std::stringstream date;
    date << std::put_time(std::localtime(&now_c), "%Y-%m-%d %X");

    file << group << "," << circuit << "," << date.str();
    for (const auto& [key, value]: times) {
        file << "," << value / 1000000.0;
    }
    file << parallel_iterations << "," << iterations << std::endl;
    file.close();
}
