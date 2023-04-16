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
#define DEBUG true
#ifdef DEBUG
#define THREAD_SAFE_PRINT(value) \
    do { \
        if (DEBUG) { \
            int tid = omp_get_thread_num(); \
            std::ofstream ofs("H:/Uni/Masterarbeit/pyzx/thesis/thread_" + std::to_string(tid) + "_output.txt", std::ios_base::app); \
            ofs << value; \
        } \
    } while(0)
#else
#define THREAD_SAFE_PRINT(value) 
#endif


namespace zx {
    

    class ExtractorParallel {
    public:
        ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, int thread_num, std::unordered_map<size_t, int>* claimed_vertices, Measurement measurement = Measurement(true));
        ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, Measurement measurement = Measurement(true));
        ExtractorParallel* other_extractor;
        std::unordered_map<size_t, int>* claimed_vertices; // Vertices marked by the thread in parallel execution
        
        void extract();
        void finalizeExtraction(std::map<zx::Qubit, zx::Vertex> other_frontier, std::unordered_set<size_t> claimed_neighbors_other);
        std::vector<size_t> frontierToInputs();
        std::map<zx::Qubit, zx::Vertex> frontier;
        bool parallelize = true;

        std::unordered_map<size_t, std::unordered_set<size_t>> deleted_edges;
        std::unordered_map<size_t, std::unordered_set<size_t>> added_edges;
        std::unordered_set<size_t> claimed_neighbors;

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

        bool extractCNOT();

        bool processFrontier();


        std::stringstream ts_printstream;

        std::vector<zx::Vertex> get_frontier_neighbors();
        std::vector<zx::Vertex> get_frontier_neighbors_parallel(std::vector<zx::Vertex>* frontier_values);

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



        void printVector(std::vector<size_t> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    void printVector(std::vector<dd::Qubit> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    void printVector(std::vector<int> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    void printVector(std::map<int, int> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v.first << " : " << v.second << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }
    void printVector(std::map<zx::Qubit, zx::Vertex> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v.first << " : " << v.second << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }
    void printVector(std::map<zx::Vertex, zx::Vertex> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v.first << " : " << v.second << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }
    void printVector(std::map<zx::Vertex, int> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v.first << " : " << v.second << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    void printVector(std::vector<bool> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( v << " ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    void printVector(std::vector<std::pair<int, int>> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( "(" << v.first << "," << v.second << ") ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    void printVector(std::set<std::pair<zx::Vertex, zx::Vertex>> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT( "(" << v.first << "," << v.second << ") ,");
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    void printMatrix(gf2Mat matrix) {
        for(auto const& row : matrix) {
            printVector(row);
        }
        THREAD_SAFE_PRINT( std::endl);
    }

    };
    
    void testParallelExtraction(std::string circuitName="vbe_adder_3.qasm", std::string measurementGroup="1", bool parallelization=false);

} // namespace zx