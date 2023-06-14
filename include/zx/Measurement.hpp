#pragma once
#include <windows.h>
#include <iostream>
#include <algorithm>
#include <cstddef>
#include <optional>
#include <map>
#include <chrono>
#include <tuple>
#include <fstream>
#include <ctime>
#include <sstream>
#include <iomanip>

class Measurement {
public:
    Measurement(bool enabled=true);
    std::map<std::string, double> times;
    bool enabled;

    void addMeasurement(std::string name, std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end);

    void printMeasurements(std::string group, std::string circuit, int parallel_iterations, int iterations, std::string filename="H:/Uni/Masterarbeit/measurements.csv");
};

