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
#include "zx/BenchmarkData.hpp"
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
#include <atomic>
#define DEBUG false
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


#define THREAD_SAFE_PRINT_2(value) \
    do { \
            int tid = omp_get_thread_num(); \
            std::ofstream ofs("H:/Uni/Masterarbeit/pyzx/thesis/thread_" + std::to_string(tid) + "_measurements.txt", std::ios_base::app); \
            ofs << value; \
    } while(0)


namespace zx {

    class ExtractorConfig {
        public:
            bool perm_optimization = true;
            bool parallel_allow_claimed_vertices_for_cnot = true;
            bool parallel_allow_cnot = true;
    };
    

    class ExtractorParallel {
    public:
        ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, int thread_num, std::unordered_map<size_t, int>* claimed_vertices, Measurement measurement = Measurement(true));
        ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, Measurement measurement = Measurement(true));

        void configure(const ExtractorConfig& config) {
            perm_optimization = config.perm_optimization;
            parallel_allow_claimed_vertices_for_cnot = config.parallel_allow_claimed_vertices_for_cnot;
            parallel_allow_cnot = config.parallel_allow_cnot;
        }

        ExtractorParallel* other_extractor;
        std::unordered_map<size_t, int>* claimed_vertices; // Vertices marked by the thread in parallel execution
        void printStatistics();
        int extract();
        int finalizeExtraction(std::map<zx::Qubit, zx::Vertex> other_frontier);
        std::vector<size_t> frontierToInputs();
        std::map<zx::Qubit, zx::Vertex> frontier;
        bool parallelize = true;

        std::unordered_map<size_t, std::unordered_set<size_t>> deleted_edges;
        std::unordered_map<size_t, std::unordered_set<size_t>> added_edges;
        std::unordered_set<size_t> claimed_neighbors;

        double time_total = 0;
        int iteration = 0;

        int parallel_iterations = 0;
        int total_iterations = 0;
        std::string circuit_name = "";
        std::string measurement_group = "";

        // Whether to use heuristic optimization for column order before gauss elim
        bool perm_optimization = false;
        bool parallel_allow_claimed_vertices_for_cnot = true;
        bool parallel_allow_cnot = true;

        // Time for extraction operations
        double time_extr_par_cnot = 0;
        double time_extr_par_cz = 0;
        double time_extr_par_fp = 0;

        double time_extr_seq_cnot = 0;
        double time_extr_seq_cz = 0;
        double time_extr_seq_fp = 0;
        
        double time_cnot_failed_extraction = 0;

        double time_cnot_gauss = 0;
        double time_cnot_biadj = 0;
        double time_cnot_optimal = 0;
        double time_cnot_neighbors = 0;

        // Number of extraction operations
        int num_extr_par_cnot = 0;
        int num_extr_par_cz = 0;
        int num_extr_par_fp = 0;

        int num_extr_seq_cnot = 0;
        int num_extr_seq_cz = 0;
        int num_extr_seq_fp = 0;

        int failedCnots = 0;

        // Number of gates created during extraction
        int num_gates_cnot = 0;
        int num_gates_cz = 0;
        int num_gates_phase = 0;
        int num_gates_h = 0;
        int num_gates_swap = 0;

        bool isFinished() const {
            return finished.load(std::memory_order::memory_order_relaxed);
        }   

        BenchmarkData ExtractorParallel::createBenchmarkData() {
            BenchmarkData data;

            // Fill the BenchmarkData object with the extractor's information
            data.time_total = time_total;
            data.parallel_iterations = parallel_iterations;
            data.total_iterations = total_iterations;
            data.time_extr_par_cnot = time_extr_par_cnot;
            data.time_extr_par_cz = time_extr_par_cz;
            data.time_extr_par_fp = time_extr_par_fp;
            data.time_extr_seq_cnot = time_extr_seq_cnot;
            data.time_extr_seq_cz = time_extr_seq_cz;
            data.time_extr_seq_fp = time_extr_seq_fp;
            data.time_cnot_failed_extraction = time_cnot_failed_extraction;
            data.time_cnot_gauss = time_cnot_gauss;
            data.time_cnot_biadj = time_cnot_biadj;
            data.time_cnot_optimal = time_cnot_optimal;
            data.time_cnot_neighbors = time_cnot_neighbors;
            data.num_extr_par_cnot = num_extr_par_cnot;
            data.num_extr_par_cz = num_extr_par_cz;
            data.num_extr_par_fp = num_extr_par_fp;
            data.num_extr_seq_cnot = num_extr_seq_cnot;
            data.num_extr_seq_cz = num_extr_seq_cz;
            data.num_extr_seq_fp = num_extr_seq_fp;
            data.failedCnots = failedCnots;
            data.num_gates_cnot = num_gates_cnot;
            data.num_gates_cz = num_gates_cz;
            data.num_gates_phase = num_gates_phase;
            data.num_gates_h = num_gates_h;
            data.num_gates_swap = num_gates_swap;

            return data;
        }
        

    private:
        qc::QuantumComputation& circuit;
        ZXDiagram& diag;
        Measurement measurement;
        std::vector<size_t> inputs;
        std::vector<size_t> outputs;
        int thread_num;
        std::atomic<bool> finished{false};
        
        
        std::vector<size_t> frontier_neighbors;
        
        void initFrontier();

        void extractRZ_CZ();

        int extractCNOT(); // 0 = success, 1 = interrupted, 2 = re-do

        bool processFrontier();

        void extractOutputHadamards();


        std::stringstream ts_printstream;

        bool get_frontier_neighbors();
        bool get_frontier_neighbors_parallel_new();
        bool get_frontier_neighbors_parallel(std::vector<zx::Vertex>* frontier_values);

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
            THREAD_SAFE_PRINT_2( v << " ,");
        }
        THREAD_SAFE_PRINT_2( std::endl);
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
            THREAD_SAFE_PRINT_2( v.first << " : " << v.second << " ,");
        }
        THREAD_SAFE_PRINT_2( std::endl);
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
            THREAD_SAFE_PRINT_2( v << " ,");
        }
        THREAD_SAFE_PRINT_2( std::endl);
    }

    void printVector(std::vector<std::pair<int, int>> vec) {
        for(auto const& v : vec) {
            THREAD_SAFE_PRINT_2( "(" << v.first << "," << v.second << ") ,");
        }
        THREAD_SAFE_PRINT_2( std::endl);
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
        THREAD_SAFE_PRINT_2( std::endl);
    }

    };
    
    BenchmarkData testParallelExtraction(std::string circuitName="vbe_adder_3.qasm", std::string measurementGroup="1", bool parallelization=false, const ExtractorConfig& config = ExtractorConfig(), bool random=false, int randomQubits=0);

} // namespace zx