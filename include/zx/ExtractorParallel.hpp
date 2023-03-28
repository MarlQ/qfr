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

    class ExtractorParallel {
    public:
        ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, int thread_num, std::map<size_t, int>* claimed_vertices, Measurement measurement = Measurement(true));

        ExtractorParallel* other_extractor;
        std::map<size_t, int>* claimed_vertices; // Vertices marked by the thread in parallel execution
        
        void extract();
        void finalizeExtraction(std::map<zx::Qubit, zx::Vertex> other_frontier);
        std::vector<size_t> frontierToInputs();
        std::map<zx::Qubit, zx::Vertex> frontier;
    private:
        qc::QuantumComputation& circuit;
        ZXDiagram& diag;
        Measurement measurement;
        std::vector<size_t> inputs;
        std::vector<size_t> outputs;
        int thread_num;
        
        
        std::vector<size_t> frontier_neighbors;
        
        void initFrontier();

        void extractRZ_CZ();

        void extractCNOT();

        void processFrontier();

        std::vector<zx::Vertex> get_frontier_neighbors();
        std::vector<zx::Vertex> get_frontier_neighbors_parallel();

        gf2Mat getAdjacencyMatrix(const std::vector<zx::Vertex>& vertices_from, const std::vector<Vertex>& vertices_to);

        std::unordered_map<int, int> column_optimal_swap(zx::gf2Mat& matrix);

        std::optional<std::unordered_map<int, int>> find_targets(std::unordered_map<int, std::unordered_set<int>> conn, std::unordered_map<int, std::unordered_set<int>> connr, std::unordered_map<int, int> target = {});

        void row_add(zx::gf2Mat& matrix, int r0, int r1, std::vector<std::pair<zx::Qubit, zx::Qubit>>& rowOperations);

        std::vector<std::pair<zx::Qubit, zx::Qubit>> gaussElimination(zx::gf2Mat& matrix);

        bool contains(std::vector<size_t> vec, zx::Vertex val);

        bool contains(std::map<zx::Qubit, zx::Vertex> vec, zx::Vertex val);

        std::vector<std::pair<int, int>>permutation_as_swaps(std::map<int, int> perm);




        void claimOutputs();
        bool claim(size_t vertex);
        bool isClaimed(size_t vertex);
        bool isClaimedBySelf(size_t vertex);
        bool isClaimedAnother(size_t vertex);





    };
    
    void testParallelExtraction(std::string circuitName="vbe_adder_3.qasm", std::string measurementGroup="1");

} // namespace zx