#define NOMINMAX
#include "zx/ExtractorParallel.hpp"

#include "algorithms/RandomCliffordCircuit.hpp"
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
#include <optional>
#include <tuple>
#define _USE_MATH_DEFINES
#include <cmath>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#define DEBUG false


/**
 * TODO-List:
 * - [x] Initialize both extractors
 * - [x] Claim inputs/outputs
 * - [x] Claim frontier
 * - [X] CZ/Phase gate extraction
 * - [ ] Check whether CNOT is necessary
 * 
 * 
 * 
 * ------- After parallel extraction
 * - [ ] 
 */

namespace zx {

    std::vector<Vertex> mapToVector(std::map<zx::Qubit, zx::Vertex>& vertices) {
        std::vector<Vertex> retVerts;

        for( std::map<zx::Qubit, zx::Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it ) {
            retVerts.push_back( it->second );
         }
        return retVerts;
    }


    

    ExtractorParallel::ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, int thread_num, std::unordered_map<size_t, int>* claimed_vertices, Measurement measurement):
        circuit(circuit), diag(diag), thread_num(thread_num), claimed_vertices(claimed_vertices), measurement(measurement), inputs(diag.getInputs()), outputs(diag.getOutputs()) {
        claimOutputs();
        omp_init_lock(&lock);
    };


    ExtractorParallel::ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, Measurement measurement):
     circuit(circuit), diag(diag), thread_num(0), parallelize(false), measurement(measurement), inputs(diag.getInputs()), outputs(diag.getOutputs()) {
        
    };

    ExtractorParallel* other_extractor;

    int ExtractorParallel::finalizeExtraction(std::map<zx::Qubit, zx::Vertex> other_frontier) {
        parallelize = false;

        // Remove all edges between the other extractor's frontier
        // and vertices that were claimed by that extractor, 
        // as those are edges that have already been handled.
        //std::cout << "Basic preparation for finalization... " << thread_num << std::endl;
        //diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test1.json", mapToVector(frontier));
                    //other_extractor->diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test2.json", mapToVector(other_extractor->frontier));
        THREAD_SAFE_PRINT( "Basic preparation for finalization..." << std::endl);
        for(auto [qubit, v] : other_frontier) {
            auto neighbors = diag.getNeighborVertices(v);

            for(auto w : neighbors) {
                if(isClaimedAnother(w)) {
                    //std::cout << "Removing edge " << "(" << v << ", " << w << ")" << std::endl;
                    diag.removeEdge(v, w);
                }
                THREAD_SAFE_PRINT( "Removing edge " << v << " , " << w << std::endl);
            }

            // Then add a simple connection to the appropriate input vertex
            auto input_vert = inputs[qubit];
            THREAD_SAFE_PRINT( "Adding edge " << v << " , " << input_vert << std::endl);
            if(!diag.connected(v, input_vert)) {
                diag.addEdge(v, input_vert);
                //std::cout << "Adding edge " << "(" << v << ", " << input_vert << ")" << std::endl; 
            } 
            else {
                diag.setEdgeType(v, input_vert, EdgeType::Simple);
                //std::cout << "Setting edge to simple " << "(" << v << ", " << input_vert << ")" << std::endl; 
            }
            // Remove phase
            if (!diag.phase(v).isZero()) {
                //std::cout << "Removing phase at " << v << std::endl; 
                THREAD_SAFE_PRINT( "Removing phase at " << v << std::endl);
                diag.setPhase(v, PiExpression());
            } 
            // Add entry to frontier if it does not exist
            if(frontier.count(qubit) == 0) {
                //std::cout << "Adding qubit " << qubit << std::endl; 
                frontier[qubit] = v;
            }     
        }
        THREAD_SAFE_PRINT( "CNOT-edges preparation for finalization..." << std::endl);

        for (const auto& entry : other_extractor->deleted_edges) {
            size_t v = entry.first;
            const std::unordered_set<size_t>& deleted_neighbors = entry.second;

            for (size_t w : deleted_neighbors) {
                // Delete the edge between v and w
                THREAD_SAFE_PRINT( "Removing edge " << v << " , " << w << std::endl);
                diag.removeEdge(v,w);
            }
        }
        for (const auto& entry : other_extractor->added_edges) {
            size_t v = entry.first;
            const std::unordered_set<size_t>& added_neighbors = entry.second;

            for (size_t w : added_neighbors) {
                // Add edge between v and w
                
                if(!diag.connected(v, w)) {
                    THREAD_SAFE_PRINT( "Adding edge " << v << " , " << w << std::endl);
                    diag.addEdge(v,w,EdgeType::Hadamard);
                }
                else {
                    THREAD_SAFE_PRINT( "Already connected " << v << " , " << w << std::endl);
                }
            }
        }
        extractOutputHadamards();

        THREAD_SAFE_PRINT( "Preparation finished." << std::endl);

        THREAD_SAFE_PRINT( "Inputs:" << std::endl);
        if(DEBUG)printVector(inputs);

        THREAD_SAFE_PRINT( "Outputs:" << std::endl);
        if(DEBUG)printVector(outputs);
        size_t originalVerts = diag.getNVertices(); 
        size_t previousVerts = originalVerts;
        int i = 0;
        while (frontier.size() > 0) {
            auto a = omp_get_wtime();
            extractRZ_CZ();
            auto b = omp_get_wtime();
            time_extr_seq_cz += (b - a) * 1000.0;

            a = omp_get_wtime();
            extractCNOT();
            b = omp_get_wtime();
            time_extr_seq_cnot += (b - a) * 1000.0;
            a = omp_get_wtime();
            processFrontier();
            b = omp_get_wtime();
            time_extr_seq_fp += (b - a) * 1000.0;

            if(DEBUG) THREAD_SAFE_PRINT( "Iteration " << i << " thread " << omp_get_thread_num() << std::endl);
            i++;

            /* if(i % 10 == 0) {
                size_t currentVerts = diag.getNVertices(); 
                float percentage = 100 - (((float) currentVerts) / ((float) originalVerts)) * 100;
                std::cout << "Completion (B): " << percentage << "%" << std::endl;
                 if(currentVerts <= previousVerts) {
                    std::cout << "STUCK " << currentVerts << " | " << previousVerts << std::endl;
                    diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test1.json", mapToVector(frontier));
                    other_extractor->diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test2.json", mapToVector(other_extractor->frontier));
                    exit(0);
                }
                previousVerts = currentVerts; *
            }*/
        };

        THREAD_SAFE_PRINT( "Finished extraction. Reversing circuit and finding swaps..." << std::endl);

        THREAD_SAFE_PRINT( "Inputs:" << std::endl);
        if(DEBUG)printVector(inputs);

        THREAD_SAFE_PRINT( "Frontier:" << std::endl);
        if(DEBUG)printVector(frontier);

        THREAD_SAFE_PRINT( "Outputs:" << std::endl);
        if(DEBUG)printVector(outputs);

        THREAD_SAFE_PRINT( "EDGES: " << std::endl);
        if (DEBUG) {
            for (auto v: diag.getEdges()) {
                std::cout << "(" << v.first << ", " << v.second << "), ";
            }
        }
        THREAD_SAFE_PRINT( std::endl);

        // Find swap gates
        std::map<int, int> swaps;
        bool               leftover_swaps = false;

        for (size_t q = 0; q < outputs.size(); ++q) {
            auto incident = diag.incidentEdges(outputs[q]);

            for (auto& edge: incident) {
                zx::Vertex w  = edge.to;
                auto       it = std::find(inputs.begin(), inputs.end(), w);
                if (it != inputs.end()) {
                    // Check if edge to input is hadamard
                    if (edge.type == zx::EdgeType::Hadamard) {
                        circuit.h(q);
                        num_gates_h++;
                        diag.setEdgeType(outputs[q], edge.to, EdgeType::Simple);
                    }
                    auto qw = (zx::Qubit)(it - inputs.begin());
                    if (((size_t)qw) != q) {
                        THREAD_SAFE_PRINT( "Found swap at " << q << " , " << qw << std::endl);
                        leftover_swaps = true;
                    }
                    swaps[q] = qw;
                    break;
                }
            }
        }
        if (leftover_swaps) {
            THREAD_SAFE_PRINT( "Creating swaps... " << std::endl);
            // Check for swaps
            for (auto s: permutation_as_swaps(swaps)) {
                circuit.swap(s.first, s.second);
                num_gates_swap++;
            }
        }

        // Reverse circuit
        circuit.reverse();

        // Additional CNOT reduction.
        ////auto //begin = std::chrono::steady_clock::now();
        qc::CircuitOptimizer::cancelCNOTs(circuit);
        //auto end = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:cancelCNOTs", begin, end);
        parallelize = true;
        return i;
    }

    void ExtractorParallel::extractOutputHadamards() {
        for (auto v: frontier) { // IMPROVE: Iterate over outputs
            auto incident = diag.incidentEdges(v.second);
            for (auto edge: incident) {
                zx::Vertex w = edge.to;

                if (!contains(outputs, w)) {
                    continue;
                }

                if (edge.type == zx::EdgeType::Hadamard) {
                    THREAD_SAFE_PRINT( "Adding Hadamard at " << v.first << std::endl);
                    circuit.h(v.first);
                    num_gates_h++;
                    diag.setEdgeType(v.second, w, zx::EdgeType::Simple);
                }
            }
        }
    }

    int ExtractorParallel::extract() {
        THREAD_SAFE_PRINT( "Extractor " << omp_get_thread_num() << " started!" << std::endl);
        if(parallelize && parallel_frontier_processing) {
            omp_set_lock(&lock);
        }
        initFrontier();
        size_t originalVerts = diag.getNVertices(); 

        extractOutputHadamards();

        iteration = 1; // Iteration counter. For debugging only.
        if(parallelize && parallel_frontier_processing) {
            omp_unset_lock(&lock);
        }
        if(parallelize && frontier.size() <= 0) {
            finished.store(true, std::memory_order::memory_order_relaxed);
            return iteration;
        }
        while (frontier.size() > 0) {
            ////auto //begin = std::chrono::steady_clock::now();
            if(parallelize && parallel_frontier_processing) {
                omp_set_lock(&lock);
            }
            auto a = omp_get_wtime();
            extractRZ_CZ();
            auto b = omp_get_wtime();
            if(parallelize) time_extr_par_cz += (b - a) * 1000.0;
            else time_extr_seq_cz += (b - a) * 1000.0;
            //auto end = std::chrono::steady_clock::now();
            //measurement.addMeasurement("extract:extractRZ_CZ", begin, end);

            //extractRZ_CZ();
/* 
            auto verts = diag.getVertices();
        std::vector<zx::Vertex> temp;

        // Check for remaining vertices
        for (auto v: frontier) {
            auto neighbors = diag.getNeighborVertices(v.second); //FIXME: Not necessary?
            THREAD_SAFE_PRINT( "neighbors of " << v.second << " :" << std::endl);
            if(DEBUG)printVector(neighbors);
        } */ 

            
            

            //begin = std::chrono::steady_clock::now();
            a = omp_get_wtime();
            int interrupted_cnot = 2;
            double X = 0.0;
            
            int u = 0;
            while(interrupted_cnot == 2) { // 0 = success, 1 = interrupted, 2 = re-do
                if(u == 1) X = omp_get_wtime();
                interrupted_cnot = extractCNOT();
                u++;
            }
            b = omp_get_wtime();
            if(parallelize) time_extr_par_cnot += (b - a) * 1000.0;
            else time_extr_seq_cnot += (b - a) * 1000.0;
            if(X > 0.0) time_cnot_failed_extraction += (b - X) * 1000.0;

            if(parallelize && parallel_frontier_processing) {
                omp_unset_lock(&lock);
            }

            //end = std::chrono::steady_clock::now();
            //measurement.addMeasurement("extract:extractCNOT", begin, end);
            //begin = std::chrono::steady_clock::now();
            if(parallelize && parallel_frontier_processing) {
                omp_set_lock(&lock);
            }
            a = omp_get_wtime();
            bool interrupted_processing = !processFrontier();
            b = omp_get_wtime();
            if(parallelize) time_extr_par_fp += (b - a) * 1000.0;
            else time_extr_seq_fp += (b - a) * 1000.0;
            //end = std::chrono::steady_clock::now();
            //measurement.addMeasurement("extract:processFrontier", begin, end);
            if(DEBUG) THREAD_SAFE_PRINT( "Iteration " << iteration << " thread " << omp_get_thread_num() << std::endl);
            
            iteration++;
/*             if(iteration % 10 == 0) {
                size_t currentVerts = diag.getNVertices(); 
                float percentage = 100 - (((float) currentVerts) / ((float) originalVerts)) * 100;
                std::cout << "Completion: " << percentage << "%" << std::endl;
            } */

            if(parallelize && parallel_frontier_processing) {
                omp_unset_lock(&lock);
            }

            if( parallelize && (interrupted_cnot == 1 ||  interrupted_processing || other_extractor->isFinished())) {
                a = omp_get_wtime();
                extractRZ_CZ();
                b = omp_get_wtime();
                if(parallelize) time_extr_par_cz += (b - a) * 1000.0;
                else time_extr_seq_cz += (b - a) * 1000.0;
                if(DEBUG) THREAD_SAFE_PRINT( "Thread " << thread_num << " was interrupted! ------------------------------------------------------------------------" << std::endl);

                finished.store(true, std::memory_order::memory_order_relaxed);
                //std::cout << "K THX BYE " << interrupted_cnot << std::endl;
                return iteration;
            }
            
        };
        if(parallelize) {
            std::cout << "THIS SHOULD NEVER HAPPEN " << std::endl;
            exit(0);  
        }

        THREAD_SAFE_PRINT( "Finished extraction. Reversing circuit and finding swaps..." << std::endl);

        THREAD_SAFE_PRINT( "Inputs:" << std::endl);
        if(DEBUG)printVector(inputs);

        THREAD_SAFE_PRINT( "Frontier:" << std::endl);
        if(DEBUG)printVector(frontier);

        THREAD_SAFE_PRINT( "Outputs:" << std::endl);
        if(DEBUG)printVector(outputs);

        THREAD_SAFE_PRINT( "EDGES: " << std::endl);
        if (DEBUG) {
            for (auto v: diag.getEdges()) {
                std::cout << "(" << v.first << ", " << v.second << "), ";
            }
        }
        THREAD_SAFE_PRINT( std::endl);

        // Find swap gates
        std::map<int, int> swaps;
        bool               leftover_swaps = false;

        for (size_t q = 0; q < outputs.size(); ++q) {
            auto incident = diag.incidentEdges(outputs[q]);

            for (auto& edge: incident) {
                zx::Vertex w  = edge.to;
                auto       it = std::find(inputs.begin(), inputs.end(), w);
                if (it != inputs.end()) {
                    // Check if edge to input is hadamard
                    if (edge.type == zx::EdgeType::Hadamard) {
                        circuit.h(q);
                        num_gates_h++;
                        diag.setEdgeType(outputs[q], edge.to, EdgeType::Simple);
                    }
                    auto qw = (zx::Qubit)(it - inputs.begin());
                    if (((size_t)qw) != q) {
                        THREAD_SAFE_PRINT( "Found swap at " << q << " , " << qw << std::endl);
                        leftover_swaps = true;
                    }
                    swaps[q] = qw;
                    break;
                }
            }
        }
        if (leftover_swaps) {
            THREAD_SAFE_PRINT( "Creating swaps... " << std::endl);
            // Check for swaps
            for (auto s: permutation_as_swaps(swaps)) {
                circuit.swap(s.first, s.second);
                num_gates_swap++;
            }
        }

        // Reverse circuit
        if (thread_num == 0) circuit.reverse();

        // Additional CNOT reduction.
        ////auto //begin = std::chrono::steady_clock::now();
        qc::CircuitOptimizer::cancelCNOTs(circuit);
        //auto end = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:cancelCNOTs", begin, end);
        return iteration;
    }

    void ExtractorParallel::printStatistics() {
        // Print the statistics
        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << "// Thread " << thread_num << std::endl;

        std::cout << "// Time for extraction operations" << std::endl;
        std::cout << "time_extr_par_cnot = " << time_extr_par_cnot  << std::endl;
        std::cout << "time_extr_par_cz = " << time_extr_par_cz  << std::endl;
        std::cout << "time_extr_par_fp = " << time_extr_par_fp  << std::endl;
        std::cout << std::endl;

        if(thread_num == 0) {
            std::cout << "time_extr_seq_cnot = " << time_extr_seq_cnot  << std::endl;
            std::cout << "time_extr_seq_cz = " << time_extr_seq_cz  << std::endl;
            std::cout << "time_extr_seq_fp = " << time_extr_seq_fp  << std::endl;
            std::cout << std::endl;
        }

        if(parallelize) {
            std::cout << "time_cnot_failed_extraction = " << time_cnot_failed_extraction  << std::endl;
            std::cout << std::endl;
        }

        std::cout << "// Number of extraction operations" << std::endl;
        std::cout << "num_extr_par_cnot = " << num_extr_par_cnot  << std::endl;
        std::cout << "num_extr_par_cz = " << num_extr_par_cz  << std::endl;
        std::cout << "num_extr_par_fp = " << num_extr_par_fp  << std::endl;
        std::cout << std::endl;

        if(thread_num == 0) {
            std::cout << "num_extr_seq_cnot = " << num_extr_seq_cnot  << std::endl;
            std::cout << "num_extr_seq_cz = " << num_extr_seq_cz  << std::endl;
            std::cout << "num_extr_seq_fp = " << num_extr_seq_fp  << std::endl;
            std::cout << std::endl;
        }

        if(parallelize) {
            std::cout << "failedCnots = " << failedCnots  << std::endl;
            std::cout << std::endl;
        }

        std::cout << "// Number of gates created during extraction" << std::endl;
        std::cout << "num_gates_cnot = " << num_gates_cnot  << std::endl;
        std::cout << "num_gates_cz = " << num_gates_cz  << std::endl;
        std::cout << "num_gates_phase = " << num_gates_phase  << std::endl;
        std::cout << "num_gates_h = " << num_gates_h  << std::endl;
        std::cout << "num_gates_swap = " << num_gates_swap  << std::endl;


        // Write the output to a CSV file
        std::ostringstream output;

        output << measurement_group;
        output << "," << circuit_name;
        output << "," << parallel_iterations;
        output << "," << total_iterations;
        output << "," << (total_iterations > 0 ?  ((double)parallel_iterations / (double)total_iterations) : 0);
        output << "," << time_extr_par_cnot;
        output << "," << time_extr_par_cz;
        output << "," << time_extr_par_fp;

        output << "," << time_cnot_failed_extraction;

        output << "," << time_extr_seq_cnot;
        output << "," << time_extr_seq_cz;
        output << "," << time_extr_seq_fp;

        output << "," << time_cnot_neighbors;
        output << "," << time_cnot_biadj;
        output << "," << time_cnot_optimal;
        output << "," << time_cnot_gauss;
        
        // Amount of extraction operations

        output << "," << num_extr_par_cnot;
        output << "," << num_extr_par_cz;
        output << "," << num_extr_par_fp;

        output << "," << failedCnots;

        output << "," << num_extr_seq_cnot;
        output << "," << num_extr_seq_cz;
        output << "," << num_extr_seq_fp;

        // Single operation times
        output << "," << (num_extr_par_cnot > 0 ? time_extr_par_cnot / num_extr_par_cnot : 0);
        output << "," << (num_extr_par_cz > 0 ? time_extr_par_cz / num_extr_par_cz : 0);
        output << "," << (num_extr_par_fp > 0 ? time_extr_par_fp/ num_extr_par_fp : 0);

        output << "," << (num_extr_seq_cnot > 0 ? time_extr_seq_cnot / num_extr_seq_cnot : 0);
        output << "," << (num_extr_seq_cz > 0 ? time_extr_seq_cz / num_extr_seq_cz : 0);
        output << "," << (num_extr_seq_fp > 0 ? time_extr_seq_fp/ num_extr_seq_fp : 0);

        output << "," << num_gates_cnot;
        output << "," << num_gates_cz;
        output << "," << num_gates_phase;
        output << "," << num_gates_h;
        output << "," << num_gates_swap;
        output << std::endl;
        
        std::ofstream csvFile("H:/Uni/Masterarbeit/statistics_complete.csv", std::ios::out | std::ios::app);
        csvFile << output.str(); 
        csvFile.close();
    }

    void ExtractorParallel::initFrontier() {
        THREAD_SAFE_PRINT( "Initializing frontier" << std::endl);
        for (size_t i = 0; i < outputs.size(); ++i) {
            auto v = diag.getNeighborVertices(outputs[i])[0];
            if (!contains(inputs, v) && claim(v)) {
                frontier[i] = v;
            }
        }
    }

    void ExtractorParallel::extractRZ_CZ() {
        THREAD_SAFE_PRINT( "Extracting RZ and CZ gates..." << std::endl);

        if(parallelize) num_extr_par_cz++;
        else num_extr_seq_cz++;
        ////auto //begin = std::chrono::steady_clock::now();
        // Extract RZ: Add phase-gate at v with phase p
        for (auto v: frontier) {
            if (!diag.phase(v.second).isZero()) {
                THREAD_SAFE_PRINT( "Adding phase gate at " << v.first << " with phase " << diag.phase(v.second).getConst().toDouble() << std::endl);
                if (!diag.phase(v.second).isConstant()) {
                    THREAD_SAFE_PRINT( "Error: phase is not constant!" << std::endl);
                    exit(0);
                }
                circuit.phase(v.first, diag.phase(v.second).getConst().toDouble());
                num_gates_phase++;
                diag.setPhase(v.second, PiExpression());
            }
        }
        //auto end = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:extractRZ_CZ:PhaseGates", begin, end);

        //begin = std::chrono::steady_clock::now();
        std::vector<std::pair<zx::Vertex, zx::Vertex>> edges_to_remove;
        std::set<zx::Vertex> handled_vertices;
        // Extract CZ
        for (auto v: frontier) { // IMPROVE: Can this be within the same for loop as the prior?
            for (zx::Edge e: diag.incidentEdges(v.second)) {
                auto w  = e.to;
                auto it = std::find_if(frontier.begin(), frontier.end(), [w](const auto& p) { return p.second == w; });
                if (it != frontier.end()) {
                    dd::Qubit qw = it->first;
                    
                    // Remove edge between v and w
                    if(handled_vertices.count(w) == 0) {
                        THREAD_SAFE_PRINT( "Adding CZ gate at " << v.first << "/" << it->first << std::endl);
                        circuit.z(v.first, dd::Control{qw});
                        num_gates_cz++;
                        edges_to_remove.push_back({v.second, w});
                    }   
                }
            }
            handled_vertices.emplace(v.second);
        }

        // Remove edge between v and w
        for (const auto& edge : edges_to_remove) {
            diag.removeEdge(edge.first, edge.second);
        }
        //end = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:extractRZ_CZ:CZGates", begin, end);
    }

    int ExtractorParallel::extractCNOT() {

        if(parallelize && !parallel_allow_cnot) return 1;
        THREAD_SAFE_PRINT( "Is CNOT extraction necessary?" << std::endl);

        std::vector<zx::Vertex> frontier_values;
        for (const auto& [key, value]: frontier) {
            frontier_values.push_back(value);
        }

        // Get frontier neighbors and check at the same time
        ////auto //begin = std::chrono::steady_clock::now();
        auto a = omp_get_wtime();
        bool cnot_necessary = parallelize ? 
            (parallel_allow_claimed_vertices_for_cnot ? get_frontier_neighbors_parallel(&frontier_values) : get_frontier_neighbors_parallel_new()) 
            : get_frontier_neighbors();
        if(!cnot_necessary) {
            THREAD_SAFE_PRINT( "No need for CNOT extraction. " << std::endl);
            return 0;
        }
        auto b = omp_get_wtime();
        time_cnot_neighbors += (b - a) * 1000.0;

        //auto end = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:extractCNOT:check", begin, end);
        THREAD_SAFE_PRINT( "Frontier CNOT extraction... " << std::endl);

        THREAD_SAFE_PRINT( "Frontier:" << std::endl);
        if(DEBUG)printVector(frontier);

        //begin = std::chrono::steady_clock::now();
        

        //frontier_neighbors = parallelize ? get_frontier_neighbors_parallel(&frontier_values) : get_frontier_neighbors();
        if(frontier_neighbors.size() <= 0) return 1;

        THREAD_SAFE_PRINT( "Frontier neighbors:" << std::endl);
        if(DEBUG)printVector(frontier_neighbors);

        // Get biadjacency matrix of frontier/neighbors
        a = omp_get_wtime();
        auto adjMatrix = getAdjacencyMatrix(frontier_values, frontier_neighbors);
        b = omp_get_wtime();
        time_cnot_biadj += (b - a) * 1000.0;
        //end            = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:extractCNOT:biadjacencyMatrix", begin, end);

        THREAD_SAFE_PRINT( "Adjacency Matrix:" << std::endl);
        if(DEBUG)printMatrix(adjMatrix);
        a = omp_get_wtime();
        if (perm_optimization) { // TODO: Measure time
            THREAD_SAFE_PRINT( "Finding optimal column swaps" << std::endl);
            std::unordered_map<int, int> perm = column_optimal_swap(adjMatrix);
            std::unordered_map<int, int> perm_swapped;
            for (const auto& [k, v]: perm) {
                perm_swapped[v] = k;
            }
            std::vector<size_t> neighbors2;
            for (size_t i = 0; i < frontier_neighbors.size(); ++i) {
                neighbors2.emplace_back(frontier_neighbors[perm_swapped[i]]);
            }
            THREAD_SAFE_PRINT( "New neighbors:" << std::endl);
            if(DEBUG)printVector(neighbors2);
            adjMatrix = getAdjacencyMatrix(frontier_values, neighbors2);
            THREAD_SAFE_PRINT( "New Adjacency Matrix:" << std::endl);
            //if(DEBUG)printMatrix(adjMatrix);
            
        }
        b = omp_get_wtime();
        time_cnot_optimal += (b - a) * 1000.0;
        // Gauss reduction on biadjacency matrix
        //begin                                                      = std::chrono::steady_clock::now();

        a = omp_get_wtime();
        std::vector<std::pair<zx::Qubit, zx::Qubit>> rowOperations = gaussElimination(adjMatrix);
        b = omp_get_wtime();
        time_cnot_gauss += (b - a) * 1000.0;
        
        
        if(parallelize) num_extr_par_cnot++;
        else num_extr_seq_cnot++;
        //end = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:extractCNOT:gaussElimination", begin, end);
        THREAD_SAFE_PRINT( "After Gauss Elimination:" << std::endl);
        if(DEBUG)printMatrix(adjMatrix);
        THREAD_SAFE_PRINT( "Row Operations:" << std::endl);
        if(DEBUG)printVector(rowOperations);

        std::vector<zx::Vertex> ws;
        bool singleOneRowExists = false;
        for (size_t i = 0; i < adjMatrix.size(); ++i) {
            int sum = 0; //, nonZero = 0;
            for (size_t j = 0; j < adjMatrix[i].size(); ++j) {
                sum += adjMatrix[i][j];
                /* if(adjMatrix[i][j]) {
                    nonZero = j;
                } */
            }
            if (sum == 1) {
                singleOneRowExists = true;
                break;
            }
        }
        //std::cout << "Vector ws:" << std::endl);
        if(DEBUG)printVector(ws);

        //begin = std::chrono::steady_clock::now();
        if (!singleOneRowExists) {
            THREAD_SAFE_PRINT( "Ws is 0" << std::endl);


            /* if(!parallelize) exit(0);
            
        
            return false;
            exit(0); */

            if(parallelize) return 1; // TODO: Parallelization not yet supported
            bool yzFound = false;
            //std::cout << "WS IS 0" << std::endl;
            //exit(0);

            //printVector(frontier);
            //printVector(other_extractor->frontier);
            //printMatrix(adjMatrix);
            //printVector(rowOperations);
            //diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test1.json", mapToVector(frontier));
            //(other_extractor->diag).toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test2.json", mapToVector(other_extractor->frontier));
            //exit(0);

            // Extract YZ-spiders
            for (auto v: frontier_neighbors) {
                auto v_neighbors = diag.getNeighborVertices(v);

                THREAD_SAFE_PRINT( "Neighbors of :" << v << std::endl);
                if(DEBUG)printVector(v_neighbors);

                //v_neighbors.erase(std::remove_if(v_neighbors.begin(), v_neighbors.end(), [diag](int x) { return diag.isInput(x); }), v_neighbors.end());

                for (auto w: v_neighbors) { // neighbors of neighbors
                    if (contains(inputs, w) || contains(outputs, w)) continue;
                    auto w_neighbors = diag.getNeighborVertices(w);

                    THREAD_SAFE_PRINT( "...Neighbors of :" << w << std::endl);
                    if(DEBUG)printVector(w_neighbors);

                    if (w_neighbors.size() == 1) { // Has no other neighbor
                        
                        THREAD_SAFE_PRINT( "Vertex with only one neighbor found: " << w << " with phase " << diag.phase(w) << std::endl);
                        THREAD_SAFE_PRINT( "This is a neighbor of frontier neighbor vert: " << v << std::endl);
                        if (contains(outputs, w)) {
                            THREAD_SAFE_PRINT( "ERROR: vertex is input!" << std::endl);
                        }

                        if (diag.phase(v).isZero()) { // Phase-gadget found
                                                      //FIXME: Or PI?
                                                      // IMPROVE: Maybe check earlier?
                            THREAD_SAFE_PRINT( "YZ Found " << v << std::endl);
                            yzFound = true;
                            size_t    corresponding_frontier_vertex;
                            zx::Qubit q;
                            bool vertex_found = false;

                            THREAD_SAFE_PRINT( "Frontier:" << std::endl);
                            if(DEBUG)printVector(frontier);
                            THREAD_SAFE_PRINT( "Searching coresponding frontier vert... "<< std::endl);
                            THREAD_SAFE_PRINT( "Neighborsss of :" << v << std::endl);
                            if(DEBUG)printVector(v_neighbors);
                            for (auto z: v_neighbors) { // Get the frontier vertex corresponding to v
                                THREAD_SAFE_PRINT( "Current " << z << std::endl);
                                auto it = std::find_if(frontier.begin(), frontier.end(), [z](const auto& p) { return p.second == z; });
                                if (it != frontier.end()) {
                                    corresponding_frontier_vertex = z;
                                    q                 = it->first;
                                    vertex_found = true;
                                    break;
                                }
                            }
                            THREAD_SAFE_PRINT( "Frontier vert is " << corresponding_frontier_vertex << std::endl);
                            if (vertex_found) { // Add hadmard - green - hadamard to output ---> green becomes frontier afterwards, then manually resolve hadamard to output
                                
                                THREAD_SAFE_PRINT( "Sanity 235"  << " == " << diag.isDeleted(235) << std::endl);
                                if(diag.isDeleted(v)) {
                                    THREAD_SAFE_PRINT( "Was DELETED "  << v << std::endl);
                                    //diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test2.json", mapToVector(frontier));
                                    exit(0);
                                }

                                if(diag.isDeleted(corresponding_frontier_vertex)) {
                                    THREAD_SAFE_PRINT( "Was DELETED "  << corresponding_frontier_vertex << std::endl);
                                    //diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test2.json", mapToVector(frontier));
                                    exit(0);
                                }

                                THREAD_SAFE_PRINT( "Performing pivot " << v << " and " << corresponding_frontier_vertex << std::endl);
                                zx::pivot(diag, v, corresponding_frontier_vertex);
                                THREAD_SAFE_PRINT( "Sanity 235"  << " == " << diag.isDeleted(235) << std::endl);

                                THREAD_SAFE_PRINT( "Deleted " << v << " == " << diag.isDeleted(v) << std::endl);
                                THREAD_SAFE_PRINT( "Deleted " << frontier[q] << " == " << diag.isDeleted(frontier[q]) << std::endl);

                                // Remove YZ-spider
                                THREAD_SAFE_PRINT( "Old spider " << frontier[q] << " AHH " << corresponding_frontier_vertex << std::endl);
                                if (frontier[q] != v) {
                                    THREAD_SAFE_PRINT( "Old spider " << frontier[q] << " != " << v << std::endl);
                                    //exit(0);
                                }
                                // frontier[q] = w;

                                //THREAD_SAFE_PRINT( "Removing YZ-spider " << v << std::endl);
                                

                                THREAD_SAFE_PRINT( "Frontier:" << std::endl);
                                if(DEBUG)printVector(frontier);

                                std::vector<size_t> test;
                                for (size_t i = 0; i < outputs.size(); ++i) {
                                    auto u = diag.getNeighborVertices(outputs[i]);
                                    THREAD_SAFE_PRINT( "Output neighbors of " << i << std::endl);
                                    if(DEBUG)printVector(u);
                                }
                                frontier[q] = diag.getNeighborVertices(outputs[q])[0]; // FIXME:
                                if(diag.isDeleted(w)) {
                                    THREAD_SAFE_PRINT( "...but it was deleted!?" << std::endl);
                                }

                                THREAD_SAFE_PRINT( "New frontier spider is " << frontier[q] << " on Qubit" << q << std::endl);

                                THREAD_SAFE_PRINT( "New frontier?:" << std::endl);
                                if(DEBUG)printVector(frontier);

                                THREAD_SAFE_PRINT( "Frontier neighbors:" << std::endl);
                                if(DEBUG)printVector(frontier_neighbors);


                            }
                            else {
                                THREAD_SAFE_PRINT( "NO FRONTIER NEIGHBOR?" << std::endl);
                            }
                        }
                    }
                }
            }
            //if (omp_get_thread_num() == 0) diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test2.json", mapToVector(frontier));
            if(yzFound) {
                THREAD_SAFE_PRINT( "All YZ-Spiders eliminated " << std::endl);
                extractOutputHadamards();
                extractRZ_CZ(); // FIXME: Hadamards to outputs?
                return extractCNOT(); // FIXME:
            }
            else {
                THREAD_SAFE_PRINT( "No YZ-Spider found " << std::endl);
                return 1;
            }
        }
        //end = std::chrono::steady_clock::now();

        //measurement.addMeasurement("extract:extractCNOT:YZSpiders", begin, end);
        //begin = std::chrono::steady_clock::now();

        // Remap row operations to vertex indices
        std::vector<std::pair<std::pair<zx::Qubit, zx::Vertex>, std::pair<zx::Qubit, zx::Vertex>>> rowOperations_indices;

        for (auto r: rowOperations) {
            /* auto r_1 = r.first;
            auto r_2 = r.second;

            if(parallelize) { // Since we previously may have removed values from the frontier, we have to find the true entries
                // IMPROVE: CHeck whether this is necessary first
                auto frontier_value_1 = frontier_values[r_1];
                auto frontier_value_2 = frontier_values[r_2];

                for (auto it = frontier.begin(); it != frontier.end(); ++it) {
                    if (it->second == frontier_value_1) r_1 = it->first;
                    else if(it->second == frontier_value_2) r_2 = it->first;
                }
            } */


            auto it            = frontier.begin();
            auto control_entry = std::next(it, r.second);
            auto target_entry  = std::next(it, r.first);

            /* auto control_qubit = control_entry->first;
            auto target_qubit  = target_entry->first;
            if (DEBUG) std::cout << " Entries " << control_qubit << "|" << control_entry->second << " , " << target_qubit << "|" << target_entry->second << std::endl;

            circuit.x(target_qubit, dd::Control{(dd::Qubit)control_qubit});
            if (DEBUG) std::cout << "Added CNOT (T:" << target_qubit << ", C:" << control_qubit << ")" << std::endl;

            auto ftarg = frontier.at(control_qubit);
            auto fcont = frontier.at(target_qubit); */



            rowOperations_indices.emplace_back(std::pair<std::pair<zx::Qubit, zx::Vertex>, std::pair<zx::Qubit, zx::Vertex>>{*target_entry, *control_entry});
        }

        if(parallelize) {
            bool conflict = false;
            std::unordered_map<size_t, std::unordered_set<size_t>> changes_deleted_edges;
            std::unordered_map<size_t, std::unordered_set<size_t>> changes_added_edges;

            //auto a = omp_get_wtime();
            #pragma omp critical(cnot) 
            {
                for (auto r: rowOperations_indices) {
                    // Check whether there are existing changes to the vertices within rowOperations

                    conflict = other_extractor->added_edges.count(r.first.second) 
                    || other_extractor->added_edges.count(r.second.second)
                    || other_extractor->deleted_edges.count(r.first.second)
                    || other_extractor->deleted_edges.count(r.second.second);

                    if(conflict) {
                        if(DEBUG) {
                            bool conflict_first = other_extractor->added_edges.count(r.first.second) 
                            || other_extractor->deleted_edges.count(r.first.second);

                            bool conflict_second = other_extractor->added_edges.count(r.second.second)
                            || other_extractor->deleted_edges.count(r.second.second);

                            if(conflict_first) THREAD_SAFE_PRINT( "Conflict found at vertex (" << r.first.second << "). Implementing changes." << std::endl);
                            if(conflict_second) THREAD_SAFE_PRINT( "Conflict found at vertex (" << r.second.second << "). Implementing changes." << std::endl);
                        }
                        
                        break;
                    }
                }

                if(conflict) {
                    // IMPROVE: Copy, implement later
                    // Implement changes
                    for(auto e : other_extractor->added_edges) {
                        for(auto to : e.second) {
                            THREAD_SAFE_PRINT( "Added edge " << e.first << " , " << to << std::endl);
                            diag.addEdge(e.first, to, EdgeType::Hadamard);
                        }

                    }
                    for(auto e : other_extractor->deleted_edges) {
                        for(auto to : e.second) {
                            THREAD_SAFE_PRINT( "Removed edge " << e.first << " , " << to << std::endl);
                            diag.removeEdge(e.first, to);
                        }

                    }

                    // Delete changes from array
                    other_extractor->added_edges.clear();
                    other_extractor->deleted_edges.clear();
                }
                else {
                    // IMPROVE: No: Add own changes to buffer array, so that the other thread is not blocked
                    // though... the performance improvement is questionable
                    THREAD_SAFE_PRINT( "No conflict!" << std::endl);
                    for (auto r: rowOperations_indices) {

                        // From row operation: r.second = r.second + r.first

                        auto it            = frontier.begin();
                        auto control_entry = r.second;
                        auto target_entry  = r.first;

                        auto control_qubit = control_entry.first;
                        auto target_qubit  = target_entry.first;
                        THREAD_SAFE_PRINT( " Entries " << control_qubit << "|" << control_entry.second << " , " << target_qubit << "|" << target_entry.second << std::endl);

                        circuit.x(target_qubit, dd::Control{(dd::Qubit)control_qubit});
                        num_gates_cnot++;
                        THREAD_SAFE_PRINT( "Added CNOT (T:" << target_qubit << ", C:" << control_qubit << ")" << std::endl);

                        auto ftarg = control_entry.second;
                        auto fcont = target_entry.second;

                        // Update diag based on row operation
                        for (auto v: diag.getNeighborVertices(fcont)) {
                            if (contains(outputs, v)) {
                                continue;
                            }

                            if (diag.connected(ftarg, v)) {
                                diag.removeEdge(ftarg, v);

                                // Remove conflicting entries
                                auto it = added_edges.find(v);
                                int removed = 0;
                                if (it != added_edges.end()) {
                                    removed = it->second.erase(ftarg); // Remove ftarg from the unordered_set
                                    if (it->second.empty())
                                        added_edges.erase(it); // Remove the entry if the unordered_set becomes empty
                                    
                                }
                                if(removed == 0) {  // DEL + ADD = NOTHING (Only remove when it was not added already)
                                    deleted_edges[v].insert(ftarg); 
                                    THREAD_SAFE_PRINT( "Removed edge (" << ftarg << ", " << v << ")" << std::endl);
                                }
                                else {
                                    THREAD_SAFE_PRINT( "Edge cancelled (" << ftarg << ", " << v << ")" << std::endl);
                                }
                                
                            } else {
                                if (contains(inputs, v)) { // v is an input
                                // FIXME: This case can not occur?
                                    THREAD_SAFE_PRINT( "Trying to remove edge but v is input " << ftarg << " vs " << v << std::endl);
                                    auto new_v = diag.insertIdentity(fcont, target_qubit, v);
                                    if (new_v) {
                                        diag.addEdge(ftarg, *new_v, zx::EdgeType::Hadamard);
                                        THREAD_SAFE_PRINT( "Added edge (" << ftarg << ", " << *new_v << ")" << std::endl);
                                    }
                                } else {
                                    diag.addEdge(ftarg, v, zx::EdgeType::Hadamard);

                                    // Remove conflicting entries
                                    auto it = deleted_edges.find(v);
                                    int removed = 0;
                                    if (it != deleted_edges.end()) {
                                        removed = it->second.erase(ftarg); // Remove ftarg from the unordered_set
                                        if (it->second.empty())
                                            deleted_edges.erase(it); // Remove the entry if the unordered_set becomes empty
                                    }
                                    if(removed == 0) { // DEL + ADD = NOTHING (Only add when it was not removed already)
                                        added_edges[v].insert(ftarg);
                                        THREAD_SAFE_PRINT( "Added edge (" << ftarg << ", " << v << ")" << std::endl);
                                    }
                                    else {
                                        THREAD_SAFE_PRINT( "Edge cancelled (" << ftarg << ", " << v << ")" << std::endl);
                                    }
                                }
                            }
                        }
                    }
                }
            } // END OF CRITICAL REGION
            //auto c = omp_get_wtime();
            //parallel_time += (c - a) * 1000.0;
            
            if(conflict) {
                THREAD_SAFE_PRINT( "Re-doing CNOT " << std::endl);
                failedCnots++;
                return 2;
            } else return 0;
        }
        

        if(!parallelize) { // TODO: Duplicate code...
            // Extract CNOTs
            for (auto r: rowOperations_indices) {
                /* auto r_1 = r.first;
                auto r_2 = r.second;

                if(parallelize) { // Since we previously may have removed values from the frontier, we have to find the true entries
                    // IMPROVE: CHeck whether this is necessary first
                    auto frontier_value_1 = frontier_values[r_1];
                    auto frontier_value_2 = frontier_values[r_2];

                    for (auto it = frontier.begin(); it != frontier.end(); ++it) {
                        if (it->second == frontier_value_1) r_1 = it->first;
                        else if(it->second == frontier_value_2) r_2 = it->first;
                    }
                    
                } */


                // From row operation: r.second = r.second + r.first

                auto it            = frontier.begin();
                auto control_entry = r.second;
                auto target_entry  = r.first;

                auto control_qubit = control_entry.first;
                auto target_qubit  = target_entry.first;
                THREAD_SAFE_PRINT( " Entries " << control_qubit << "|" << control_entry.second << " , " << target_qubit << "|" << target_entry.second << std::endl);

                circuit.x(target_qubit, dd::Control{(dd::Qubit)control_qubit});
                num_gates_cnot++;
                THREAD_SAFE_PRINT( "Added CNOT (T:" << target_qubit << ", C:" << control_qubit << ")" << std::endl);

                auto ftarg = control_entry.second;
                auto fcont = target_entry.second;

                // Update diag based on row operation
                for (auto v: diag.getNeighborVertices(fcont)) {
                    if (contains(outputs, v)) {
                        continue;
                    }

                    if (diag.connected(ftarg, v)) {
                        diag.removeEdge(ftarg, v);
                        THREAD_SAFE_PRINT( "Removed edge (" << ftarg << ", " << v << ")" << std::endl);
                    } else {
                        if (contains(inputs, v)) { // v is an input
                            THREAD_SAFE_PRINT( "Trying to remove edge but v is input " << ftarg << " vs " << v << std::endl);
                            auto new_v = diag.insertIdentity(fcont, target_qubit, v);
                            if (new_v) {
                                diag.addEdge(ftarg, *new_v, zx::EdgeType::Hadamard);
                                THREAD_SAFE_PRINT( "Added edge (" << ftarg << ", " << *new_v << ")" << std::endl);
                            }
                        } else {
                            diag.addEdge(ftarg, v, zx::EdgeType::Hadamard);
                            THREAD_SAFE_PRINT( "Added edge (" << ftarg << ", " << v << ")" << std::endl);
                        }
                    }
                }
            }
        }
        //end = std::chrono::steady_clock::now();
        //measurement.addMeasurement("extract:extractCNOT:CNOTFromOperations", begin, end);

        return 0;
    }

    bool ExtractorParallel::processFrontier() {
        THREAD_SAFE_PRINT( "Processing Frontier... " << std::endl);
        THREAD_SAFE_PRINT( "Frontier:" << std::endl);
        //if(DEBUG)printVector(frontier);

        if(parallelize) num_extr_par_fp++;
        else num_extr_seq_fp++;

        std::map<zx::Qubit, int> new_frontier;

        for (auto const& v: frontier) {
            THREAD_SAFE_PRINT( "Vertex: " << v.second << std::endl);
            auto current_neighbors = diag.getNeighborVertices(v.second);
            THREAD_SAFE_PRINT( "Neighb.:" << std::endl);
            if(DEBUG)printVector(current_neighbors);
            if (current_neighbors.size() > 2 || contains(inputs, v.second)) {
                continue; // Process later
            }
            zx::Vertex output;

            // Eliminate Hadamard chains
            bool                    uneven_hadamard = false;
            zx::Vertex              previous_vertex;
            zx::Vertex              current_vertex = v.second;
            std::vector<zx::Vertex> chain;

            bool claimed = false;

            for (zx::Vertex n: current_neighbors) {
                if (contains(outputs, n)) {
                    previous_vertex = n;
                    output          = n;
                    break;
                }
            }

            while (true) {
                zx::Vertex next_vertex;
                bool found_neighbor = false;

                for (auto n: current_neighbors) {
                    if (n != previous_vertex) {
                        if(parallelize) {
                            if(parallel_frontier_processing) {
                                if(!claim(n)) { // Take the other thread's frontier vertex
                                    if(thread_num == 1) continue; // TODO: Allow thread 1 to do this, too...
                                    omp_set_lock(&(other_extractor->lock));
                                    int qb = other_extractor->removeFromFrontier(n);
                                    
                                    if(qb >= 0) { // Take it
                                        forceClaim(n);
                                        auto neighbors = diag.getNeighborVertices(n);

                                        for(auto w : neighbors) {
                                            if(isClaimedAnother(w)) {
                                                diag.removeEdge(n, w);
                                            }
                                        }

                                        // Then add a simple connection to the appropriate input vertex
                                        auto input_vert = inputs[qb];
                                        if(!diag.connected(n, input_vert)) {
                                            diag.addEdge(n, input_vert);    
                                        }


                                        // Remove phase
                                        if (!diag.phase(n).isZero()) {
                                            diag.setPhase(n, PiExpression());
                                        } 

                                        // TODO: Only implement changes relevant to the particular vertex
                                        for(auto e : other_extractor->deleted_edges) {
                                            for(auto to : e.second) {
                                                //if(to == n) {
                                                    THREAD_SAFE_PRINT( "Removed edge " << e.first << " , " << to << std::endl);
                                                    diag.removeEdge(e.first, to);
                                                    //std::cout << "RRemoved (" << e.first << " " << to << ")" << std::endl;
                                                //}

                                            }

                                        }

                                        for(auto e : other_extractor->added_edges) {
                                            for(auto to : e.second) {
                                                //if(to == n) {
                                                    THREAD_SAFE_PRINT( "Added edge " << e.first << " , " << to << std::endl);
                                                    if(!diag.connected(e.first, to)) {
                                                        diag.addEdge(e.first, to, EdgeType::Hadamard);
                                                        //std::cout << "AAdded (" << e.first << " " << to << ")" << std::endl; 
                                                    }
                                                //}
                                            }
                                        }

                                        // Delete changes from array
                                        other_extractor->added_edges.clear();
                                        other_extractor->deleted_edges.clear();

                                        omp_unset_lock(&(other_extractor->lock)); // I believe I can't do this earlier, because then the other thread might think the vertex is up for grabs
                                    }
                                    else {
                                        omp_unset_lock(&(other_extractor->lock));
                                        continue; // Was not a frontier vertex
                                    }
                                }
                                
                            }
                            else {
                                if(!claim(n)) continue;
                            }
                        }
                        
                        next_vertex = n;
                        found_neighbor = true;
                        break;
                    }
                }
                if(!found_neighbor) {
                    THREAD_SAFE_PRINT( "No unclaimed neighbor: Stopping! " << std::endl);
                    break;
                }

                chain.emplace_back(next_vertex);

                auto edge = diag.getEdge(current_vertex, next_vertex);

                if (!edge) return false;

                if (edge->type == EdgeType::Hadamard) {
                    uneven_hadamard = !uneven_hadamard;
                }
                if (!diag.phase(next_vertex).isZero()) { // IMPROVE: Should this really break?
                    break;
                }
                if (diag.getNeighborVertices(next_vertex).size() > 2 || contains(inputs, next_vertex)) {
                    break;
                }

                previous_vertex   = current_vertex;
                current_vertex    = next_vertex;
                current_neighbors = diag.getNeighborVertices(current_vertex);
            } // END: while

            if(chain.size() <= 0) {
                THREAD_SAFE_PRINT( "No chain found." << std::endl);
                continue;
            }

            if (uneven_hadamard) {
                circuit.h(v.first);
                num_gates_h++;
                THREAD_SAFE_PRINT( "Adding Hadamard at " << v.first << std::endl);
            }
            THREAD_SAFE_PRINT( "Chain found:" << std::endl);
            if(DEBUG)printVector(chain);
            
            for (unsigned int i = 0; i < chain.size() - 1; i++) {
                THREAD_SAFE_PRINT( "Removing Vertex: " << chain[i] << std::endl);
                diag.removeVertex(chain[i]);
                if(parallelize/*  && parallel_allow_cnot && parallel_allow_claimed_vertices_for_cnot */) { // Vertex will be behind the frontier. We no longer need to store edge information about it.
                    THREAD_SAFE_PRINT( "Deleting entry " << chain[i] << std::endl);
                    #pragma omp critical(cnot) 
                    {
                        deleted_edges.erase(chain[i]);
                        added_edges.erase(chain[i]);
                    }
                }
            }
            auto edgeType      = diag.getEdge(v.second, output)->type; //CHECK: isn't this always Simple, as we removed hadamards at the start?
            auto last_in_chain = chain[chain.size() - 1];
            //auto v_qubit = v.first;
            THREAD_SAFE_PRINT( "Edge type is " << ((edgeType == EdgeType::Simple) ? "S" : "H") << std::endl);
            THREAD_SAFE_PRINT( "Removing Vertex: " << v.second << std::endl);
            diag.removeVertex(v.second);
            THREAD_SAFE_PRINT( "Adding Edge: (" << last_in_chain << "," << output << ")" << std::endl);
            diag.addEdge(last_in_chain, output, edgeType);

            if(parallelize /* && parallel_allow_cnot && parallel_allow_claimed_vertices_for_cnot */) { // Vertex will be in the frontier. We no longer need to store edge information about it.
                THREAD_SAFE_PRINT( "Deleting entry " << last_in_chain << std::endl);
                #pragma omp critical(cnot) 
                {
                    deleted_edges.erase(last_in_chain);
                    added_edges.erase(last_in_chain);
                }
            }

            if (!contains(inputs, last_in_chain)) {
                new_frontier[v.first] = (int)last_in_chain;
            } else {
                new_frontier[v.first] = -1;
            }
            //std::cout << "Frontier Changes:" << std::endl);
            //if(DEBUG)printVector(new_frontier);
        }
        THREAD_SAFE_PRINT( "Old Frontier:" << std::endl);
        if(DEBUG)printVector(frontier);

        if (new_frontier.size() > 0) {
            THREAD_SAFE_PRINT( "Frontier Changes:" << std::endl);
            if(DEBUG)printVector(new_frontier);
            for (auto entry: new_frontier) {
                //frontier.erase(std::remove(frontier.begin(), frontier.end(), entry.first), frontier.end());
                frontier.erase(entry.first);
                if (entry.second != -1) {
                    frontier[entry.first] = (zx::Vertex)entry.second;
                }
            }
        } else return false;
        THREAD_SAFE_PRINT( "New Frontier:" << std::endl);
        if(DEBUG)printVector(frontier);
        return true;
    }

    // Check whether all frontier vertices have more than two neighbors
    // and at the same time, adds non-output neighbors to frontier_neighbors
    bool ExtractorParallel::get_frontier_neighbors() {
        std::vector<Vertex> neighbors;

        for (auto v: frontier) {
            auto v_neighbors = diag.getNeighborVertices(v.second);
            if(v_neighbors.size() <= 2) return false; // CNOT not necessary: vertex only has one non-boundary neighbor
            for (auto w: v_neighbors) {
                if (!contains(outputs, w) && !contains(neighbors, w)) {
                    neighbors.emplace_back(w);
                }
            }
        }

        frontier_neighbors = neighbors;
        return true;
    }


    bool ExtractorParallel::get_frontier_neighbors_parallel_new() {
        std::vector<Vertex> neighbors;
        #pragma omp critical(claim)
        {
            for (auto v: frontier) {
                auto v_neighbors = diag.getNeighborVertices(v.second);
                int count = 0;

                for (auto w: v_neighbors) {
                    if (!contains(outputs, w) && !contains(neighbors, w) /* && ! contains(frontier, w) */) {
                        // Check whether neighbor is unmarked
                        auto it_c = claimed_vertices->find(w);
                        bool claimed_by_another = it_c != claimed_vertices->end() && it_c->second != thread_num;
                        neighbors.emplace_back(w);
                        if(!claimed_by_another) { 
                            count++;
                        }
                    }
                }


                // Check whether CNOT is necessary
                if(count <= 2) {
                    THREAD_SAFE_PRINT("Vertex with only one neighbor found: " << v.second << std::endl);
                    if(DEBUG) printVector(v_neighbors);
                    return false; // CNOT not necessary: vertex only has one non-boundary neighbor 
                }
            }
        }
        frontier_neighbors = neighbors;
        return true;
    }


    bool ExtractorParallel::get_frontier_neighbors_parallel(std::vector<zx::Vertex>* frontier_values) {
        std::vector<Vertex> neighbors;
        #pragma omp critical(claim)
        {
            for (auto it = frontier_values->cbegin(); it != frontier_values->cend();) {
                auto v_neighbors = diag.getNeighborVertices(*it);
                //std::vector<Vertex> neighbors_current;
                bool claimed_vertex_found = false;
                int neighbors_count = 0;

                for (auto w: v_neighbors) {
                    if (!contains(outputs, w) /* && ! contains(frontier, w) */) {

                        // Check whether neighbor is unclaimed
                        auto it_c = claimed_vertices->find(w);
                        bool claimed_by_another = it_c != claimed_vertices->end() && it_c->second != thread_num;
                        if(claimed_by_another) { 
                            THREAD_SAFE_PRINT("Found claimed neighbor: " << w << std::endl);
                            claimed_vertex_found = true;   
                        }
                        
                        neighbors_count++;

                        if(!contains(neighbors, w)) {
                            neighbors.emplace_back(w);
                        } 
                    }
                }

                // Check whether CNOT is necessary
                if(neighbors_count == 1 && !claimed_vertex_found) {
                    THREAD_SAFE_PRINT("Vertex with only one neighbor found: " << *it << std::endl);
                    if(DEBUG) printVector(v_neighbors);
                    return false; // CNOT not necessary: vertex only has one non-boundary neighbor 
                }
                ++it;
            }
        }
        frontier_neighbors = neighbors;
        return true;
    }

    gf2Mat ExtractorParallel::getAdjacencyMatrix(const std::vector<zx::Vertex>& vertices_from, const std::vector<Vertex>& vertices_to) {
        THREAD_SAFE_PRINT( " Matrix: " << vertices_from.size() << " x " << vertices_to.size() << std::endl);
        gf2Mat adjMat{vertices_from.size(), gf2Vec(vertices_to.size(), false)};
        for (size_t i = 0; i < vertices_from.size(); ++i) {
            for (size_t j = 0; j < vertices_to.size(); ++j) {
                if (diag.connected(vertices_from[i], vertices_to[j])) {
                    adjMat[i][j] = true;
                } else {
                    adjMat[i][j] = false;
                }
            }
        }
        return adjMat;
    }

    /** Column optimal swap function is a depth-first search method from PyZX, to reduce the number
     * of row operations, and consequently the number of CNOT gates generated.
     * It is very, very costly in terms of computation.
    */
    std::unordered_map<int, int> ExtractorParallel::column_optimal_swap(zx::gf2Mat& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();

        std::unordered_map<int, std::unordered_set<int>> connections, connectionsr;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (matrix[i][j]) {
                    connections[i].insert(j);
                    connectionsr[j].insert(i);
                }
            }
        }

        // connections: row -> column indices with non-zero elements
        // connections: column -> row indices with non-zero elements
        std::optional<std::unordered_map<int, int>> target_opt = find_targets(connections, connectionsr);
        std::unordered_map<int, int>                target;
        if (!target_opt)
            target = std::unordered_map<int, int>();
        else
            target = *target_opt;

        std::unordered_set<int> target_values;
        for (const auto& elem: target) {
            target_values.insert(elem.second);
        }

        std::vector<int> left;
        for (int i = 0; i < cols; i++) {
            if (target_values.find(i) == target_values.end()) {
                left.push_back(i);
            }
        }

        std::vector<int> right;
        for (int i = 0; i < cols; i++) {
            if (target.find(i) == target.end()) {
                right.push_back(i);
            }
        }

        auto left_it = left.begin(), right_it = right.begin();
        while (left_it != left.end() && right_it != right.end()) {
            target[*right_it] = *left_it;
            ++left_it;
            ++right_it;
        }

        return target;
    }

    std::optional<std::unordered_map<int, int>> ExtractorParallel::find_targets(std::unordered_map<int, std::unordered_set<int>> conn, std::unordered_map<int, std::unordered_set<int>> connr, std::unordered_map<int, int> target) {
        int r = conn.size();
        //int c = connr.size();

        std::unordered_set<int> claimedcols;
        std::unordered_set<int> claimedrows;

        for (auto [key, value]: target) {
            claimedcols.insert(key);
            claimedrows.insert(value);
        }

        while (true) {
            int                     min_index = -1;
            bool                    did_break = false;
            std::unordered_set<int> min_options;

            for (int i = 0; i < 1000; ++i) {
                min_options.insert(i);
            }

            for (int i = 0; i < r; ++i) {
                if (claimedrows.count(i)) continue;

                std::unordered_set<int> s; // The free columns
                for (auto j: conn[i]) {
                    if (!claimedcols.count(j)) s.insert(j);
                }

                if (s.size() == 1) {
                    int j     = *s.begin();
                    target[j] = i;
                    claimedcols.insert(j);
                    claimedrows.insert(i);
                    did_break = true;
                    break;
                }

                if (s.empty()) return std::nullopt; // contradiction

                bool found_col = false;
                for (auto j: s) {
                    std::unordered_set<int> t;
                    for (auto k: connr[j]) {
                        if (!claimedrows.count(k)) t.insert(k);
                    }
                    if (t.size() == 1) { // j can only be connected to i
                        target[j] = i;
                        claimedcols.insert(j);
                        claimedrows.insert(i);
                        found_col = true;
                        break;
                    }
                }
                if (found_col) {
                    did_break = true;
                    break;
                }
                if (s.size() < min_options.size()) {
                    min_index   = i;
                    min_options = s;
                }
            }
            if (!did_break) {                            // Didn't find any forced choices
                if (conn.size() == claimedrows.size()) { // Equivalent to if not (conn.keys() - claimedrows): ?
                    return target;                       // we are done
                }
                // XXX: --------------
                std::unordered_set<int> conn_values;
                for (auto [key, value]: conn) {
                    conn_values.insert(key);
                }

                if (conn_values == claimedrows) { // we are done
                    std::cout << "THERE IS A BUG" << std::endl;
                    exit(0);
                    return target;
                }
                // ------------------

                if (min_index == -1) {
                    std::cout << "This shouldn't happen ever" << std::endl;
                    exit(0);
                }

                // Start depth-first search
                std::unordered_map<int, int> tgt = target;
                for (int k: min_options) {
                    tgt[k]          = min_index;
                    auto new_target = find_targets(conn, connr, tgt);
                    if (new_target) return new_target;
                }
                return target;
            }
        }
    }

    void ExtractorParallel::row_add(zx::gf2Mat& matrix, int r0, int r1, std::vector<std::pair<zx::Qubit, zx::Qubit>>& rowOperations) {
        int cols = matrix[0].size();
        for (int k = 0; k < cols; ++k) {
            matrix[r1][k] = matrix[r1][k] ^ matrix[r0][k];
        }
        rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{r0, r1});
    }

    std::vector<std::pair<zx::Qubit, zx::Qubit>> ExtractorParallel::gaussElimination(zx::gf2Mat& matrix) {
        int                                          rows      = matrix.size();
        int                                          cols      = matrix[0].size();
        int                                          blocksize = 6;
        int                                          pivot_row = 0;
        std::vector<std::pair<zx::Qubit, zx::Qubit>> rowOperations;
        std::vector<int>                             pivot_cols;

        for (int sec = 0; sec < std::ceil(cols / (double)blocksize); ++sec) {
            int i0 = sec * blocksize;
            int i1 = std::min(cols, (sec + 1) * blocksize);

            std::unordered_map<std::vector<bool>, int> chunks;
            for (int r = pivot_row; r < rows; ++r) {
                std::vector<bool> t;
                t.reserve(i1 - i0);
                std::copy(matrix[r].begin() + i0, matrix[r].begin() + i1, std::back_inserter(t));
                if (std::none_of(t.begin(), t.end(), [](bool b) { return b; })) {
                    continue;
                }
                auto it = chunks.find(t);
                if (it != chunks.end()) {
                    row_add(matrix, chunks[t], r, rowOperations);
                } else {
                    chunks[t] = r;
                }
            }
            int p = i0;
            while (p < i1) {
                for (int r0 = pivot_row; r0 < rows; ++r0) {
                    if (matrix[r0][p]) {
                        if (r0 != pivot_row) {
                            row_add(matrix, r0, pivot_row, rowOperations);
                        }

                        for (int r1 = pivot_row + 1; r1 < rows; ++r1) {
                            if (pivot_row != r1 && matrix[r1][p]) {
                                row_add(matrix, pivot_row, r1, rowOperations);
                            }
                        }

                        pivot_cols.emplace_back(p);
                        pivot_row++;
                        break;
                    }
                }
                p++;
            }
        }

        pivot_row--;
        std::vector<int> pivot_cols1(pivot_cols);

        for (int sec = std::ceil(cols / (double)blocksize) - 1; sec >= 0; sec--) {
            int i0 = sec * blocksize;
            int i1 = std::min(cols, (sec + 1) * blocksize);

            std::unordered_map<std::vector<bool>, int> chunks;
            for (int r = pivot_row - 1; r >= 0; r--) {
                std::vector<bool> t;
                t.reserve(i1 - i0);
                std::copy(matrix[r].begin() + i0, matrix[r].begin() + i1, std::back_inserter(t));
                if (std::none_of(t.begin(), t.end(), [](bool b) { return b; })) {
                    continue;
                }
                auto it = chunks.find(t);
                if (it != chunks.end()) {
                    row_add(matrix, chunks[t], r, rowOperations);
                } else {
                    chunks[t] = r;
                }
            }

            while (pivot_cols1.size() != 0 && i0 <= pivot_cols1.back() < i1) {
                int pcol = pivot_cols1.back();
                pivot_cols1.pop_back();
                for (int r = 0; r < pivot_row; ++r) {
                    if (matrix[r][pcol]) {
                        row_add(matrix, pivot_row, r, rowOperations);
                    }
                }
                pivot_row--;
            }
        }
        return rowOperations;
    }

    bool ExtractorParallel::contains(std::vector<size_t> vec, zx::Vertex val) { // IMPROVE: isIn is already implemented in ZXDiagram
        return std::find(vec.begin(), vec.end(), val) != vec.end();
    }

    bool ExtractorParallel::contains(std::map<zx::Qubit, zx::Vertex> vec, zx::Vertex val) {
        return std::find_if(vec.begin(), vec.end(), [val](const auto& p) { return p.second == val; }) != vec.end();
    }

    std::vector<std::pair<int, int>> ExtractorParallel::permutation_as_swaps(std::map<int, int> perm) {
        THREAD_SAFE_PRINT( "Perm:" << std::endl);
        //if(DEBUG)printVector(perm);

        std::vector<std::pair<int, int>> swaps;
        THREAD_SAFE_PRINT( "Size: " << perm.size() << std::endl);
        THREAD_SAFE_PRINT( "Size: " << perm.size() << std::endl);

        std::vector<int> l;
        for (size_t i = 0; i < perm.size(); ++i) {
            l.emplace_back(perm[i]);
        }

        THREAD_SAFE_PRINT( "l:" << std::endl);
        //if(DEBUG)printVector(l);

        std::map<int, int> pinv;
        for (auto i: perm) {
            pinv[i.second] = i.first;
        }

        THREAD_SAFE_PRINT( "pinv:" << std::endl);
        //if(DEBUG)printVector(pinv);

        std::vector<int> linv;
        for (size_t i = 0; i < pinv.size(); ++i) {
            linv.emplace_back(pinv[i]);
        }

        THREAD_SAFE_PRINT( "linv:" << std::endl);
        //if(DEBUG)printVector(linv);

        for (size_t i = 0; i < perm.size(); ++i) {
            if ((size_t)l[i] == i) continue;
            int t1 = l[i];
            int t2 = linv[i];
            THREAD_SAFE_PRINT( "Adding swap gate at " << i << " , " << t2 << std::endl);
            swaps.emplace_back(std::pair<int, int>(i, t2));
            l[t2]    = t1;
            linv[t1] = t2;
        }
        return swaps;
    }

    // [CRITICAL] Checks whether a vertex has been claimed by the current thread
    bool ExtractorParallel::isClaimedBySelf(size_t vertex) {
        bool isClaimed = false;
        #pragma omp critical(claim)
        {
            auto it = claimed_vertices->find(vertex);
            if (it != claimed_vertices->end() && it->second == thread_num) {
                isClaimed = true;
            }
        }
        return isClaimed;
    }

    // [CRITICAL] Checks whether the vertex has been claimed by any thread
    bool ExtractorParallel::isClaimed(size_t vertex) {
        bool isClaimed = false;
        #pragma omp critical(claim)
        {
            isClaimed = (claimed_vertices->count(vertex) != 0);
        }
        return isClaimed;
    }

    // [CRITICAL] Tries to claim a vertex with the current thread. 
    // Returns true when the vertex was successfully claimed,
    // or false if it was already claimed.
    bool ExtractorParallel::claim(size_t vertex) {
        if(!parallelize) return true;
        bool success = true;
        #pragma omp critical(claim)
        {
            auto it = claimed_vertices->find(vertex);
            if(it != claimed_vertices->end() && it->second != thread_num) {
                THREAD_SAFE_PRINT( thread_num << ": Failed to claim vertex " << vertex << std::endl);
                success = false;
            }
            if(success) {
                claimed_vertices->emplace(vertex, thread_num);
                THREAD_SAFE_PRINT( thread_num << ": Successfully claimed vertex " << vertex << std::endl);
            }
            
        }
        return success;
    }

    void ExtractorParallel::forceClaim(size_t vertex) {
        if(!parallelize) return;
        #pragma omp critical(claim)
        {
            (*claimed_vertices)[vertex] = thread_num;
            THREAD_SAFE_PRINT( thread_num << ": Successfully claimed vertex " << vertex << std::endl);
        }
    }

    // [CRITICAL] Checks whether a vertex has been claimed by another thread
    bool ExtractorParallel::isClaimedAnother(size_t vertex) {  
        if(!parallelize) {
            auto it = claimed_vertices->find(vertex);
            if (it != claimed_vertices->end() && it->second != thread_num) {
                return true;
            }
            return false;
        }
        else {
            bool claimed = false;
            #pragma omp critical(claim)
            {
                auto it = claimed_vertices->find(vertex);
                if (it != claimed_vertices->end() && it->second != thread_num) {
                    claimed = true;
                }
            }
            return claimed;
        }
    }

    // [CRITICAL] Claims the outputs with the current thread
    void ExtractorParallel::claimOutputs() {
        THREAD_SAFE_PRINT( "Claiming outputs" << std::endl);
        for(auto v : outputs) {
            claim(v);
        }
    }

    /* std::vector<Vertex> ExtractorParallel::getClaimedNeighborVertices(const Vertex v) { 
        std::vector<Vertex> ret;
        const auto& incident = diag.incidentEdges(v);
        for(auto const& e : incident) {
            if(isClaimedBySelf(e.to)) {
                ret.emplace_back(e.to);
            }         
        }
        return ret;
    } */

    std::vector<size_t> ExtractorParallel::frontierToInputs() {
        return inputs;
        // TODO:
        size_t qubits = outputs.size();
        std::vector<size_t> convertedInputs;
        for(size_t i = 0; i < qubits; ++i) {
            auto it = frontier.find(i);

            if (it == frontier.end()) { // qubit not in frontier
                THREAD_SAFE_PRINT( "ERROR: Qubit not in frontier: " << i << std::endl);
                exit(-1);
            }
            size_t vertex = it->second;

            convertedInputs.emplace_back(vertex);
        }

        return convertedInputs;
    }

    BenchmarkData testParallelExtraction(std::string circuitName, std::string measurementGroup, bool parallelization, const ExtractorConfig& config, bool random, int randomQubits) {
        //std::cout << "Extraction testing" << std::endl;
        Measurement measurement;

        if (DEBUG) std::cout << "Setting up...\n";
        std::unique_ptr<qc::QuantumComputation> qc;
        if(random) {
            std::cout << "Generating random circuit with " << randomQubits << " qubits." << std::endl;
            const dd::QubitCount nq = randomQubits;

            auto dd = std::make_unique<dd::Package<>>(nq);
            qc = std::make_unique<qc::QuantumComputation>(qc::RandomCliffordCircuit(nq, nq * nq, 12345));
        }
        else {
            std::filesystem::path currentPath = std::filesystem::current_path();
            std::string relativePath = circuitName;
            std::filesystem::path outputPath = currentPath / relativePath;

            qc = std::make_unique<qc::QuantumComputation>();
            //std::cout << "Circuit " << circuitName << ":" << std::endl;
            qc->import(outputPath.string());
            std::cout << "Loading circuit " << outputPath << std::endl;
        }
        
        //std::cout << "Circuit loaded" << std::endl;
        //qc->dump("H:/Uni/Masterarbeit/pyzx/thesis/original.qasm");
        //std::cout << "ASd " << qc->getNops() << std::endl;
        //qc->printStatistics(std::cout);
        //if (DEBUG) std::cout << "Circuit to extract:" << std::endl;
        //if (DEBUG) std::cout << qc << std::endl;

        zx::ZXDiagram zxDiag = zx::FunctionalityConstruction::buildFunctionality(qc.get());

        zxDiag.toGraphlike();
/*         auto vertNumBefore = zxDiag.getNVertices();
        std::cout << "Vertices (before) " << zxDiag.getNVertices() << std::endl;
        std::cout << "Edges (before) " << zxDiag.getNEdges() << std::endl; */

        //std::cout << "Simplifying" << std::endl;
        zx::interiorCliffordSimp(zxDiag);
        //zx::fullReduce(zxDiag);

/* 
        auto vertNumAfter = zxDiag.getNVertices();
        std::cout << "Vertices (after) " << zxDiag.getNVertices() << std::endl;
        std::cout << "Edges (after) " << zxDiag.getNEdges() << std::endl;

        auto percentageReduction = vertNumBefore / vertNumAfter * 100;
        std::cout << "% reduction " << percentageReduction << std::endl; */

        //std::cout << "Extracting" << std::endl;
        int iterations = 0;
        int parallel_iterations = 0;
        BenchmarkData data;
        if(parallelization) {
            //std::cout << "Starting parallel extraction" << std::endl;

            //auto begin = std::chrono::steady_clock::now(); // Start measurement
            auto a = omp_get_wtime();

            qc::QuantumComputation qc_extracted = qc::QuantumComputation(zxDiag.getNQubits());
            zx::ZXDiagram          zxDiag_reversed = zxDiag.reverse();
            qc::QuantumComputation qc_extracted_2  = qc::QuantumComputation(zxDiag_reversed.getNQubits());

            std::unordered_map<size_t, int> claimed_vertices;

            ExtractorParallel extractor1(qc_extracted, zxDiag, 0, &claimed_vertices, measurement);
            ExtractorParallel extractor2(qc_extracted_2, zxDiag_reversed, 1, &claimed_vertices, measurement);
            extractor1.configure(config);
            extractor2.configure(config);
            extractor1.circuit_name = circuitName;
            extractor2.circuit_name = circuitName;
            extractor1.measurement_group = measurementGroup;
            extractor2.measurement_group = measurementGroup;

            extractor1.other_extractor = &extractor2; // IMPROVE: needed?
            extractor2.other_extractor = &extractor1;

            int iterations_1 = 0;
            int iterations_2 = 0;

            #pragma omp parallel num_threads(2) shared(claimed_vertices)
            {
                if (omp_get_thread_num() == 0) {
                    iterations_1 = extractor1.extract();
                    if(DEBUG) std::cout << "Extractor 0 finished" << std::endl;
                } else {
                    iterations_2 = extractor2.extract();
                    if(DEBUG) std::cout << "Extractor 1 finished" << std::endl;
                }
            }
            auto end_parallel = omp_get_wtime();
            extractor1.time_parallel = (end_parallel - a) * 1000.0;
            //zxDiag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test1.json", mapToVector(extractor1.frontier));
            //zxDiag_reversed.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test2.json", mapToVector(extractor1.frontier));

            // Extract rest
            //THREAD_SAFE_PRINT( "Finalizing extraction" << std::endl);
            //THREAD_SAFE_PRINT( "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl);
            //THREAD_SAFE_PRINT( "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl);
            //THREAD_SAFE_PRINT( "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl);
            /* extractor1.deleted_edges = extractor2.deleted_edges;
            extractor1.added_edges = extractor2.added_edges; */
            int iterations_seq = extractor1.finalizeExtraction(extractor2.frontier);

            // Combine diagrams
            //THREAD_SAFE_PRINT( "Combining diagrams" << std::endl);
            qc_extracted_2.combine(qc_extracted);

            //auto end = std::chrono::steady_clock::now(); // End measurement
            //measurement.addMeasurement("extract", begin, end);

            auto b = omp_get_wtime();
            extractor1.time_total = (b - a) * 1000.0;

            //qc_extracted.dump("H:/Uni/Masterarbeit/pyzx/thesis/extracted2.qasm");

            //std::cout << "Finished Circuit" << std::endl;
            //std::cout << qc_extracted_2 << std::endl;

            //qc_extracted_2.dump("H:/Uni/Masterarbeit/pyzx/thesis/extracted.qasm");
            //std::cout << "Circuit " << circuitName << std::endl;
            //std::cout << "Parallel iterations: " << iterations_1 << "||" << iterations_2 << std::endl;

            iterations = std::max(iterations_1, iterations_2) + iterations_seq;
            parallel_iterations = std::min(iterations_1, iterations_2);

            float degree_of_parallelism = std::min(iterations_1, iterations_2) / iterations;
            extractor1.iteration_diff = std::max(iterations_1, iterations_2) - std::min(iterations_1, iterations_2);
            extractor1.parallel_iterations = parallel_iterations;
            extractor1.total_iterations = iterations;
            extractor2.total_iterations = 0;
            extractor2.parallel_iterations = 0;
            //extractor1.printStatistics();
            //extractor2.printStatistics();
            //qc_extracted_2.printStatistics(std::cout);

            data = extractor1.createBenchmarkData();
            BenchmarkData data2 = extractor2.createBenchmarkData();
            data.averageData(data2);
        }
        else {
            //std::cout << "Starting non-parallel extraction" << std::endl;

            //auto begin = std::chrono::steady_clock::now(); // Start measurement
            auto a = omp_get_wtime();

            qc::QuantumComputation qc_extracted = qc::QuantumComputation(zxDiag.getNQubits());
            ExtractorParallel extractor1(qc_extracted, zxDiag, measurement);
            extractor1.configure(config);
            extractor1.circuit_name = circuitName;
            extractor1.measurement_group = measurementGroup;
            iterations = extractor1.extract();

            //auto end = std::chrono::steady_clock::now(); // End measurement
            //measurement.addMeasurement("extract", begin, end);
            auto b = omp_get_wtime();
            //double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
            extractor1.time_total = (b - a) * 1000.0;
            //extractor1.time_total = duration / 1000000.0;
            
            //std::cout << "Finished Circuit" << std::endl;
            //std::cout << qc_extracted << std::endl;

            //qc_extracted.dump("H:/Uni/Masterarbeit/pyzx/thesis/extracted.qasm");
            extractor1.total_iterations = iterations;
            //extractor1.printStatistics();
            //qc_extracted.printStatistics(std::cout);
            //std::cout << "Circuit " << circuitName << ":" << std::endl;
            data = extractor1.createBenchmarkData();
        }
        //std::cout << "Total iterations: " << iterations << std::endl;
        
        //std::cout << "Circuit to extract:" << std::endl;
        //std::cout << qc << std::endl;
       
        
        //measurement.printMeasurements(measurementGroup, circuitName, parallel_iterations, iterations, "H:/Uni/Masterarbeit/measurements_new.csv");
        return data;
    }
} // namespace zx