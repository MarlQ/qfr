#pragma once

#define NOMINMAX
#include "CircuitOptimizer.hpp"
#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "Rational.hpp"
#include "Simplify.hpp"
#include "Utils.hpp"
#include "ZXDiagram.hpp"
#include "dd/Control.hpp"
#include "zx/FunctionalityConstruction.hpp"
#include "zx/Measurement.hpp"
//#include "EquivalenceCheckingManager.hpp"
//#include "gtest/gtest.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <map>
#include <omp.h>
#include <optional>
#include <tuple>

namespace zx {

    class Extractor {
    public:
        Extractor(qc::QuantumComputation& circuit, ZXDiagram& diag, Measurement measurement = Measurement(true), bool parallelize = false);

        Extractor* other_extractor;

        void extract();

    private:
        qc::QuantumComputation& circuit;
        ZXDiagram& diag;
        Measurement measurement;
        bool parallelize;
        std::vector<size_t> inputs;
        std::vector<size_t> outputs;

        std::map<zx::Qubit, zx::Vertex> frontier;
        std::vector<size_t> frontier_neighbors;

        std::map<zx::Qubit, zx::Vertex> initFrontier();

        void extractRZ_CZ();

        void extractCNOT();

        void processFrontier();

        void removeParallelOverlap();

        std::vector<zx::Vertex> get_frontier_neighbors();

        gf2Mat getAdjacencyMatrix(const std::vector<zx::Vertex>& vertices_from, const std::vector<Vertex>& vertices_to);

        std::unordered_map<int, int> column_optimal_swap(zx::gf2Mat& matrix);

        std::optional<std::unordered_map<int, int>> find_targets(std::unordered_map<int, std::unordered_set<int>> conn, std::unordered_map<int, std::unordered_set<int>> connr, std::unordered_map<int, int> target = {});

        void row_add(zx::gf2Mat& matrix, int r0, int r1, std::vector<std::pair<zx::Qubit, zx::Qubit>>& rowOperations);

        std::vector<std::pair<zx::Qubit, zx::Qubit>> gaussElimination(zx::gf2Mat& matrix);

        bool contains(std::vector<size_t> vec, zx::Vertex val);

        bool contains(std::map<zx::Qubit, zx::Vertex> vec, zx::Vertex val);

        std::vector<std::pair<int, int>>permutation_as_swaps(std::map<int, int> perm);
    };
    
    void testExtraction(std::string circuitName="vbe_adder_3.qasm", std::string measurementGroup="1");

} // namespace zx