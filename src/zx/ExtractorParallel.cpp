#define NOMINMAX
#include "zx/ExtractorParallel.hpp"

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
#define DEBUG true

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


    void printVector(std::vector<size_t> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }

    void printVector(std::vector<dd::Qubit> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }

    void printVector(std::vector<int> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }

    void printVector(std::map<int, int> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v.first << " : " << v.second << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }
    void printVector(std::map<zx::Qubit, zx::Vertex> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v.first << " : " << v.second << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }
    void printVector(std::map<zx::Vertex, zx::Vertex> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v.first << " : " << v.second << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }
    void printVector(std::map<zx::Vertex, int> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v.first << " : " << v.second << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }

    void printVector(std::vector<bool> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << v << " ,";
        }
        if(DEBUG)std::cout << std::endl;
    }

    void printVector(std::vector<std::pair<int, int>> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << "(" << v.first << "," << v.second << ") ,";
        }
        if(DEBUG)std::cout << std::endl;
    }

    void printVector(std::set<std::pair<zx::Vertex, zx::Vertex>> vec) {
        for(auto const& v : vec) {
            if(DEBUG)std::cout << "(" << v.first << "," << v.second << ") ,";
        }
        if(DEBUG)std::cout << std::endl;
    }

    void printMatrix(gf2Mat matrix) {
        for(auto const& row : matrix) {
            printVector(row);
        }
        if(DEBUG)std::cout << std::endl;
    }

    ExtractorParallel::ExtractorParallel(qc::QuantumComputation& circuit, ZXDiagram& diag, int thread_num, std::map<size_t, int>* claimed_vertices, Measurement measurement):
        circuit(circuit), diag(diag), thread_num(thread_num), claimed_vertices(claimed_vertices), measurement(measurement), inputs(diag.getInputs()), outputs(diag.getOutputs()) {
        claimOutputs();
    };

    ExtractorParallel* other_extractor;

    void ExtractorParallel::finalizeExtraction(std::map<zx::Qubit, zx::Vertex> other_frontier) {
        

        // Remove all edges between the other extractor's frontier
        // and vertices that were claimed by that extractor, 
        // as those are edges that have already been handled.
        for(auto [qubit, v] : other_frontier) {
            auto neighbors = diag.getNeighborVertices(v);

            for(auto w : neighbors) {
                if(!isClaimedAnother(w)) continue;
                
                diag.removeEdge(v, w);
            }

            // Then add a connection to the appropriate input vertex
            auto input_vert = inputs[qubit];
            diag.addEdge(v, input_vert);     

            // Remove phase
            if (!diag.phase(v).isZero()) {
                if (DEBUG) std::cout << "Removing phase at " << v << std::endl;
                diag.setPhase(v, PiExpression());
            }      

        }

        if (DEBUG) std::cout << "Inputs:" << std::endl;
        if(DEBUG)printVector(inputs);

        if (DEBUG) std::cout << "Outputs:" << std::endl;
        if(DEBUG)printVector(outputs);

        int i = 0;
        while (frontier.size() > 0) {

            extractRZ_CZ();
            extractCNOT();
            processFrontier();
            if(DEBUG) std::cout << "Iteration " << i << " thread " << omp_get_thread_num() << std::endl;
            i++;
        };

        if (DEBUG) std::cout << "Finished extraction. Reversing circuit and finding swaps..." << std::endl;

        if (DEBUG) std::cout << "Inputs:" << std::endl;
        if(DEBUG)printVector(inputs);

        if (DEBUG) std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        if (DEBUG) std::cout << "Outputs:" << std::endl;
        if(DEBUG)printVector(outputs);

        if (DEBUG) std::cout << "EDGES: " << std::endl;
        if (DEBUG) {
            for (auto v: diag.getEdges()) {
                std::cout << "(" << v.first << ", " << v.second << "), ";
            }
        }
        if (DEBUG) std::cout << std::endl;

        // Find swap gates
        // TODO: Measure time
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
                        diag.setEdgeType(outputs[q], edge.to, EdgeType::Simple);
                    }
                    auto qw = (zx::Qubit)(it - inputs.begin());
                    if (((size_t)qw) != q) {
                        if (DEBUG) std::cout << "Found swap at " << q << " , " << qw << std::endl;
                        leftover_swaps = true;
                    }
                    swaps[q] = qw;
                    break;
                }
            }
        }
        if (leftover_swaps) {
            if (DEBUG) std::cout << "Creating swaps... " << std::endl;
            // Check for swaps
            for (auto s: permutation_as_swaps(swaps)) {
                circuit.swap(s.first, s.second);
            }
        }

        // Reverse circuit
        circuit.reverse();

        // Additional CNOT reduction.
        auto begin = std::chrono::steady_clock::now();
        qc::CircuitOptimizer::cancelCNOTs(circuit);
        auto end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:cancelCNOTs", begin, end);
        return;
    }

    void ExtractorParallel::extract() {
        std::cout << "Extractor " << omp_get_thread_num() << " started!" << std::endl;
        initFrontier();
        

        for (auto v: frontier) { // IMPROVE: Iterate over outputs
            auto incident = diag.incidentEdges(v.second);
            for (auto edge: incident) {
                zx::Vertex w = edge.to;

                if (!contains(outputs, w)) {
                    continue;
                }

                if (edge.type == zx::EdgeType::Hadamard) {
                    if (DEBUG) std::cout << "Adding Hadamard at " << v.first << std::endl;
                    circuit.h(v.first);
                    diag.setEdgeType(v.second, w, zx::EdgeType::Simple);
                }
            }
        }

        int i = 1; // Iteration counter. For debugging only.

        while (frontier.size() > 0) {

            auto begin = std::chrono::steady_clock::now();
            extractRZ_CZ();
            auto end = std::chrono::steady_clock::now();
            measurement.addMeasurement("extract:extractRZ_CZ", begin, end);



/* 
            auto verts = diag.getVertices();
        std::vector<zx::Vertex> temp;

        // Check for remaining vertices
        for (auto v: frontier) {
            auto neighbors = diag.getNeighborVertices(v.second); //FIXME: Not necessary?
            if (DEBUG) std::cout << "neighbors of " << v.second << " :" << std::endl;
            if(DEBUG)printVector(neighbors);
        } */


return;

            

            begin = std::chrono::steady_clock::now();
            extractCNOT();
            end = std::chrono::steady_clock::now();
            measurement.addMeasurement("extract:extractCNOT", begin, end);

            begin = std::chrono::steady_clock::now();
            processFrontier();
            end = std::chrono::steady_clock::now();
            measurement.addMeasurement("extract:processFrontier", begin, end);
            if(DEBUG) std::cout << "Iteration " << i << " thread " << omp_get_thread_num() << std::endl;
            i++;
        };

        if (DEBUG) std::cout << "Finished extraction. Reversing circuit and finding swaps..." << std::endl;

        if (DEBUG) std::cout << "Inputs:" << std::endl;
        if(DEBUG)printVector(inputs);

        if (DEBUG) std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        if (DEBUG) std::cout << "Outputs:" << std::endl;
        if(DEBUG)printVector(outputs);

        if (DEBUG) std::cout << "EDGES: " << std::endl;
        if (DEBUG) {
            for (auto v: diag.getEdges()) {
                std::cout << "(" << v.first << ", " << v.second << "), ";
            }
        }
        if (DEBUG) std::cout << std::endl;

        // Find swap gates
        // TODO: Measure time
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
                        diag.setEdgeType(outputs[q], edge.to, EdgeType::Simple);
                    }
                    auto qw = (zx::Qubit)(it - inputs.begin());
                    if (((size_t)qw) != q) {
                        if (DEBUG) std::cout << "Found swap at " << q << " , " << qw << std::endl;
                        leftover_swaps = true;
                    }
                    swaps[q] = qw;
                    break;
                }
            }
        }
        if (leftover_swaps) {
            if (DEBUG) std::cout << "Creating swaps... " << std::endl;
            // Check for swaps
            for (auto s: permutation_as_swaps(swaps)) {
                circuit.swap(s.first, s.second);
            }
        }

        // Reverse circuit
        if (omp_get_thread_num() == 0) circuit.reverse();

        // Additional CNOT reduction.
        auto begin = std::chrono::steady_clock::now();
        qc::CircuitOptimizer::cancelCNOTs(circuit);
        auto end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:cancelCNOTs", begin, end);
        return;
    }

    void ExtractorParallel::initFrontier() {
        std::cout << "Initializing frontier" << std::endl;
        for (size_t i = 0; i < outputs.size(); ++i) {
            auto v = diag.getNeighborVertices(outputs[i])[0];
            if (!contains(inputs, v) && claim(v)) {
                frontier[i] = v;
            }
        }
    }

    void ExtractorParallel::extractRZ_CZ() {
        if (DEBUG) std::cout << "Extracting RZ and CZ gates..." << std::endl;

        auto begin = std::chrono::steady_clock::now();
        // Extract RZ: Add phase-gate at v with phase p
        for (auto v: frontier) {
            if (!diag.phase(v.second).isZero()) {
                if (DEBUG) std::cout << "Adding phase gate at " << v.first << " with phase " << diag.phase(v.second).getConst().toDouble() << std::endl;
                if (!diag.phase(v.second).isConstant()) {
                    if (DEBUG) std::cout << "Error: phase is not constant!" << std::endl;
                    exit(0);
                }
                circuit.phase(v.first, diag.phase(v.second).getConst().toDouble());
                diag.setPhase(v.second, PiExpression());
            }
        }
        auto end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:extractRZ_CZ:PhaseGates", begin, end);

        begin = std::chrono::steady_clock::now();

        // Extract CZ
        for (auto v: frontier) { // IMPROVE: Can this be within the same for loop as the prior?
            for (zx::Edge e: diag.incidentEdges(v.second)) {
                auto w  = e.to;
                auto it = std::find_if(frontier.begin(), frontier.end(), [w](const auto& p) { return p.second == w; });
                if (it != frontier.end()) {
                    dd::Qubit qw = it->first;
                    if (DEBUG) std::cout << "Adding CZ gate at " << v.first << "/" << it->first << std::endl;

                    circuit.z(v.first, dd::Control{qw});

                    // Remove edge between v and w
                    diag.removeEdge(v.second, w);
                }
            }
        }
        end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:extractRZ_CZ:CZGates", begin, end);
    }

    void ExtractorParallel::extractCNOT() {
        if (DEBUG) std::cout << "Is CNOT extraction necessary?" << std::endl;
        //bool vertsRemaining = false;
        auto verts = diag.getVertices();

        std::vector<zx::Vertex> temp;

        // Check for remaining vertices
        auto begin = std::chrono::steady_clock::now();
        for (auto v: frontier) {
            auto neighbors = diag.getNeighborVertices(v.second); //FIXME: Not necessary?
            if (DEBUG) std::cout << "neighbors of " << v.second << " :" << std::endl;
            if(DEBUG)printVector(neighbors);
            if (neighbors.size() <= 2) {
                if (DEBUG) std::cout << "No need for CNOT extraction. " << std::endl;
                return;
            }
        }
        auto end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:extractCNOT:check", begin, end);
        if (DEBUG) std::cout << "Frontier CNOT extraction... " << std::endl;

        if (DEBUG) std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        if (DEBUG) std::cout << "Verts remaining:" << std::endl;
        if(DEBUG)printVector(temp);
        begin = std::chrono::steady_clock::now();
        frontier_neighbors = get_frontier_neighbors();

        if (DEBUG) std::cout << "Frontier neighbors:" << std::endl;
        if(DEBUG)printVector(frontier_neighbors);

        // Get biadjacency matrix of frontier/neighbors
        std::vector<zx::Vertex> frontier_values;
        for (const auto& [key, value]: frontier) {
            frontier_values.push_back(value);
        }
        // TODO: Omit frontier vertices that are neighbors to marked vertices
        auto adjMatrix = getAdjacencyMatrix(frontier_values, frontier_neighbors);
        end            = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:extractCNOT:biadjacencyMatrix", begin, end);

        if (DEBUG) std::cout << "Adjacency Matrix:" << std::endl;
        if(DEBUG)printMatrix(adjMatrix);

        bool perm_optimization = true;
        if (perm_optimization) { // TODO: Measure time
            if (DEBUG) std::cout << "Finding optimal column swaps" << std::endl;
            std::unordered_map<int, int> perm = column_optimal_swap(adjMatrix);
            std::unordered_map<int, int> perm_swapped;
            for (const auto& [k, v]: perm) {
                perm_swapped[v] = k;
            }
            std::vector<size_t> neighbors2;
            for (size_t i = 0; i < frontier_neighbors.size(); ++i) {
                neighbors2.emplace_back(frontier_neighbors[perm_swapped[i]]);
            }
            if (DEBUG) std::cout << "New neighbors:" << std::endl;
            if(DEBUG)printVector(neighbors2);
            adjMatrix = getAdjacencyMatrix(frontier_values, neighbors2);
            if (DEBUG) std::cout << "New Adjacency Matrix:" << std::endl;
            //if(DEBUG)printMatrix(adjMatrix);
        }

        // Gauss reduction on biadjacency matrix
        begin                                                      = std::chrono::steady_clock::now();
        std::vector<std::pair<zx::Qubit, zx::Qubit>> rowOperations = gaussElimination(adjMatrix);

        end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:extractCNOT:gaussElimination", begin, end);
        if (DEBUG) std::cout << "After Gauss Elimination:" << std::endl;
        if(DEBUG)printMatrix(adjMatrix);
        if (DEBUG) std::cout << "Row Operations:" << std::endl;
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
        //std::cout << "Vector ws:" << std::endl;
        if(DEBUG)printVector(ws);

        begin = std::chrono::steady_clock::now();
        if (!singleOneRowExists) {
            if (DEBUG) std::cout << "Ws is 0" << std::endl;
            diag.toJSON("H:\\Uni\\Masterarbeit\\pyzx\\thesis\\test.json", inputs);
            exit(0);

            // Extract YZ-spiders
            for (auto v: frontier_neighbors) {
                auto v_neighbors = diag.getNeighborVertices(v);

                //v_neighbors.erase(std::remove_if(v_neighbors.begin(), v_neighbors.end(), [diag](int x) { return diag.isInput(x); }), v_neighbors.end());

                for (auto w: v_neighbors) {
                    if (contains(inputs, w) || contains(outputs, w)) continue;
                    auto w_neighbors = diag.getNeighborVertices(w);

                    if (w_neighbors.size() == 1) {
                        if (DEBUG) std::cout << "Vertex with only one neighbor found: " << w << " with phase " << diag.phase(w) << std::endl;
                        if (contains(outputs, w)) {
                            if (DEBUG) std::cout << "ERROR: vertex is input!" << std::endl;
                        }

                        if (diag.phase(v).isZero()) { // Phase-gadget found
                                                      //FIXME: Or PI?

                            size_t    frontier_neighbor;
                            zx::Qubit q;

                            for (auto z: v_neighbors) {
                                auto it = std::find_if(frontier.begin(), frontier.end(), [z](const auto& p) { return p.second == z; });
                                if (it != frontier.end()) {
                                    frontier_neighbor = z;
                                    q                 = it->first;
                                    break;
                                }
                            }
                            if (frontier_neighbor) { // Add hadmard - green - hadamard to output ---> green becomes frontier afterwards, then manually resolve hadamard to output
                                zx::pivot(diag, v, frontier_neighbor);

                                // Remove YZ-spider
                                if (DEBUG) std::cout << "Old spider " << frontier[q] << std::endl;
                                if (frontier[q] != v) {
                                    if (DEBUG) std::cout << "Old spider " << frontier[q] << " != " << v << std::endl;
                                    //exit(0);
                                }
                                frontier[q] = w;

                                //if(DEBUG)std::cout << "Removing YZ-spider " << v << std::endl;
                                if (DEBUG) std::cout << "New frontier spider is " << w << " on Qubit" << q << std::endl;
                            }
                        }
                    }
                }
            }
            return;
        }
        end = std::chrono::steady_clock::now();

        // IMPROVE: Removing duplicate row operations
        //filter_duplicate_cnots(rowOperations);

        measurement.addMeasurement("extract:extractCNOT:YZSpiders", begin, end);
        begin = std::chrono::steady_clock::now();

        // Extract CNOTs
        for (auto r: rowOperations) {
            // From row operation: r.second = r.second + r.first

            auto it            = frontier.begin();
            auto control_entry = std::next(it, r.second);
            auto target_entry  = std::next(it, r.first);

            auto control_qubit = control_entry->first;
            auto target_qubit  = target_entry->first;
            if (DEBUG) std::cout << " Entries " << control_qubit << "|" << control_entry->second << " , " << target_qubit << "|" << target_entry->second << std::endl;

            circuit.x(target_qubit, dd::Control{(dd::Qubit)control_qubit});
            if (DEBUG) std::cout << "Added CNOT (T:" << target_qubit << ", C:" << control_qubit << ")" << std::endl;

            auto ftarg = frontier.at(control_qubit);
            auto fcont = frontier.at(target_qubit);

            // Update diag based on row operation
            for (auto v: diag.getNeighborVertices(fcont)) {
                if (contains(outputs, v)) {
                    continue;
                }

                if (diag.connected(ftarg, v)) {
                    diag.removeEdge(ftarg, v);
                    if (DEBUG) std::cout << "Removed edge (" << ftarg << ", " << v << ")" << std::endl;
                } else {
                    if (contains(inputs, v)) { // v is an input
                        if (DEBUG) std::cout << "Trying to remove edge but v is input " << ftarg << " vs " << v << std::endl;
                        auto new_v = diag.insertIdentity(fcont, target_qubit, v);
                        if (new_v) {
                            diag.addEdge(ftarg, *new_v, zx::EdgeType::Hadamard);
                            if (DEBUG) std::cout << "Added edge (" << ftarg << ", " << *new_v << ")" << std::endl;
                        }
                    } else {
                        diag.addEdge(ftarg, v, zx::EdgeType::Hadamard);
                        if (DEBUG) std::cout << "Added edge (" << ftarg << ", " << v << ")" << std::endl;
                    }
                }
            }
        }
        end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract:extractCNOT:CNOTFromOperations", begin, end);

        return;
    }

    void ExtractorParallel::processFrontier() {
        if (DEBUG) std::cout << "Processing Frontier... " << std::endl;
        if (DEBUG) std::cout << "Frontier:" << std::endl;
        //if(DEBUG)printVector(frontier);

        std::map<zx::Qubit, int> new_frontier;

        for (auto const& v: frontier) {
            if (DEBUG) std::cout << "Vertex: " << v.second << std::endl;
            auto current_neighbors = diag.getNeighborVertices(v.second);
            if (DEBUG) std::cout << "Neighb.:" << std::endl;
            //if(DEBUG)printVector(current_neighbors);
            if (current_neighbors.size() > 2 || contains(inputs, v.second)) {
                continue; // Process later
            }
            zx::Vertex output;

            // Eliminate Hadamard chains
            bool                    uneven_hadamard = false;
            zx::Vertex              previous_vertex;
            zx::Vertex              current_vertex = v.second;
            std::vector<zx::Vertex> chain;

            for (zx::Vertex n: current_neighbors) {
                if (contains(outputs, n)) {
                    previous_vertex = n;
                    output          = n;
                    break;
                }
            }

            while (true) {
                zx::Vertex next_vertex;

                for (auto n: current_neighbors) {
                    if (n != previous_vertex) {
                        next_vertex = n;
                    }
                }

                chain.emplace_back(next_vertex);

                auto edge = diag.getEdge(current_vertex, next_vertex);

                if (!edge) return;

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
            }

            if (uneven_hadamard) {
                circuit.h(v.first);
                if (DEBUG) std::cout << "Adding Hadamard at " << v.first << std::endl;
            }
            if (DEBUG) std::cout << "Chain found:" << std::endl;
            //if(DEBUG)printVector(chain);

            for (unsigned int i = 0; i < chain.size() - 1; i++) {
                if (DEBUG) std::cout << "Removing Vertex: " << chain[i] << std::endl;
                diag.removeVertex(chain[i]);
            }

            auto edgeType      = diag.getEdge(v.second, output)->type; //CHECK: isn't this always Simple, as we removed hadamards at the start?
            auto last_in_chain = chain[chain.size() - 1];
            //auto v_qubit = v.first;
            if (DEBUG) std::cout << "Edge type is " << ((edgeType == EdgeType::Simple) ? "S" : "H") << std::endl;
            if (DEBUG) std::cout << "Removing Vertex: " << v.second << std::endl;
            diag.removeVertex(v.second);
            if (DEBUG) std::cout << "Adding Edge: (" << last_in_chain << "," << output << ")" << std::endl;
            diag.addEdge(last_in_chain, output, edgeType);

            if (!contains(inputs, last_in_chain)) {
                new_frontier[v.first] = (int)last_in_chain;
            } else {
                new_frontier[v.first] = -1;
            }
            //std::cout << "Frontier Changes:" << std::endl;
            //if(DEBUG)printVector(new_frontier);
        }
        if (DEBUG) std::cout << "Old Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        if (new_frontier.size() > 0) {
            if (DEBUG) std::cout << "Frontier Changes:" << std::endl;
            if(DEBUG)printVector(new_frontier);
            for (auto entry: new_frontier) {
                //frontier.erase(std::remove(frontier.begin(), frontier.end(), entry.first), frontier.end());
                frontier.erase(entry.first);
                if (entry.second != -1) {
                    frontier[entry.first] = (zx::Vertex)entry.second;
                }
            }
        }
        if (DEBUG) std::cout << "New Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);
    }

    std::vector<zx::Vertex> ExtractorParallel::get_frontier_neighbors() {
        std::vector<Vertex> neighbors;

        for (auto v: frontier) {
            auto v_neighbors = diag.getNeighborVertices(v.second);
            for (auto w: v_neighbors) {
                if (!contains(outputs, w) && !contains(neighbors, w)) {
                    neighbors.emplace_back(w);
                }
            }
        }

        return neighbors;
    }


    // In contrast to the normal get_frontier_neighbors, this function only adds unmarked neighbors.
    // If a neighbor was marked already, it retroactively removes the froniter vertex, and all its neighbors
    std::vector<zx::Vertex> ExtractorParallel::get_frontier_neighbors_parallel() {
        std::vector<Vertex> neighbors;

        for (auto it = frontier.cbegin(); it != frontier.cend();) {
            auto v_neighbors = diag.getNeighborVertices(it->second);
            std::vector<Vertex> neighbors_current;
            bool marked_found = false;

            for (auto w: v_neighbors) {
                if (!contains(outputs, w) && !contains(neighbors, w) && ! contains(frontier, w)) {

                    // Check whether neighbor is unmarked
                    bool marked = claimed_vertices->count(w) != 0;
                    if(marked) { // Marked found: stop
                      marked_found = true;   
                      break;
                    }
                    
                    neighbors_current.emplace_back(w);
                }
            }

            if(marked_found) {
              // Found a single marked vertex: Remove frontier vertex
              std::cout << "Found marked neighbor, removing " << it->second << " , " << omp_get_thread_num() << std::endl;
              frontier.erase(it++);
            }
            else {
              // All unmarked: add to frontier neighbors and mark all
              for(auto w: neighbors_current) {
                std::cout << "FN Marked " << w << " , " << omp_get_thread_num() << std::endl;
                neighbors.emplace_back(w);
                claimed_vertices->emplace(w, omp_get_thread_num());
              }
              
              ++it;
            }

        }

        return neighbors;
    }

    gf2Mat ExtractorParallel::getAdjacencyMatrix(const std::vector<zx::Vertex>& vertices_from, const std::vector<Vertex>& vertices_to) {
        if (DEBUG) std::cout << " Matrix: " << vertices_from.size() << " x " << vertices_to.size() << std::endl;
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
        if (DEBUG) std::cout << "Perm:" << std::endl;
        //if(DEBUG)printVector(perm);

        std::vector<std::pair<int, int>> swaps;
        if (DEBUG) std::cout << "Size: " << perm.size() << std::endl;
        if (DEBUG) std::cout << "Size: " << perm.size() << std::endl;

        std::vector<int> l;
        for (size_t i = 0; i < perm.size(); ++i) {
            l.emplace_back(perm[i]);
        }

        if (DEBUG) std::cout << "l:" << std::endl;
        //if(DEBUG)printVector(l);

        std::map<int, int> pinv;
        for (auto i: perm) {
            pinv[i.second] = i.first;
        }

        if (DEBUG) std::cout << "pinv:" << std::endl;
        //if(DEBUG)printVector(pinv);

        std::vector<int> linv;
        for (size_t i = 0; i < pinv.size(); ++i) {
            linv.emplace_back(pinv[i]);
        }

        if (DEBUG) std::cout << "linv:" << std::endl;
        //if(DEBUG)printVector(linv);

        for (size_t i = 0; i < perm.size(); ++i) {
            if ((size_t)l[i] == i) continue;
            int t1 = l[i];
            int t2 = linv[i];
            if (DEBUG) std::cout << "Adding swap gate at " << i << " , " << t2 << std::endl;
            swaps.emplace_back(std::pair<int, int>(i, t2));
            l[t2]    = t1;
            linv[t1] = t2;
        }
        return swaps;
    }

    // [CRITICAL] Checks whether a vertex has been claimed by the current thread
    bool ExtractorParallel::isClaimedBySelf(size_t vertex) {
        #pragma omp critical(claim)
        {
            auto it = claimed_vertices->find(vertex);
            if (it != claimed_vertices->end() && it->second == thread_num) {
                return true;
            }
            return false;
        }
    }

    // [CRITICAL] Checks whether the vertex has been claimed by any thread
    bool ExtractorParallel::isClaimed(size_t vertex) {
        #pragma omp critical(claim)
        {
            return (claimed_vertices->count(vertex) != 0);
        }
    }

    // [CRITICAL] Tries to claim a vertex with the current thread. 
    // Returns true when the vertex was successfully claimed,
    // or false if it was already claimed.
    bool ExtractorParallel::claim(size_t vertex) {
        #pragma omp critical(claim)
        {
            if(claimed_vertices->count(vertex) != 0) {
                std::cout << thread_num << ": Failed to claim vertex " << vertex << std::endl;
                return false;
            }
            claimed_vertices->emplace(vertex, thread_num);
            std::cout << thread_num << ": Successfully claimed vertex " << vertex << std::endl;
            return true;
        }
    }

    // [CRITICAL] Checks whether a vertex has been claimed by another thread
    bool ExtractorParallel::isClaimedAnother(size_t vertex) {
        #pragma omp critical(claim)
        {
            auto it = claimed_vertices->find(vertex);
            if (it != claimed_vertices->end() && it->second != thread_num) {
                return true;
            }
            return false;
        }
    }

    // [CRITICAL] Claims the outputs with the current thread
    void ExtractorParallel::claimOutputs() {
        std::cout << "Claiming outputs" << std::endl;
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
                std::cout << "ERROR: Qubit not in frontier: " << i << std::endl;
                exit(-1);
            }
            size_t vertex = it->second;

            convertedInputs.emplace_back(vertex);
        }

        return convertedInputs;
    }

    

    void testParallelExtraction(std::string circuitName, std::string measurementGroup) {
        std::cout << "Parallel Extraction testing" << std::endl;
        bool        parallelize = true;
        Measurement measurement;

        if (DEBUG) std::cout << "Setting up...\n";
        qc::QuantumComputation qc{};
        std::cout << "Circuit " << circuitName << ":" << std::endl;
        qc.import("H:/Uni/Masterarbeit/qcec/" + circuitName);

        qc.dump("H:/Uni/Masterarbeit/pyzx/thesis/original.qasm");

        if (DEBUG) std::cout << "Circuit to extract:" << std::endl;
        if (DEBUG) std::cout << qc << std::endl;
        auto begin = std::chrono::steady_clock::now();

        zx::ZXDiagram zxDiag = zx::FunctionalityConstruction::buildFunctionality(&qc);
        auto          end    = std::chrono::steady_clock::now();
        measurement.addMeasurement("buildFunctionality", begin, end);

        zxDiag.toGraphlike();

        std::cout << "Simplifying" << std::endl;
        begin = std::chrono::steady_clock::now();
        zx::interiorCliffordSimp(zxDiag);
        //zx::fullReduce(zxDiag);
        end = std::chrono::steady_clock::now();
        measurement.addMeasurement("interiorCliffordSimp", begin, end);

        std::cout << "Extracting" << std::endl;

        begin = std::chrono::steady_clock::now();

        qc::QuantumComputation qc_extracted = qc::QuantumComputation(zxDiag.getNQubits());

        zx::ZXDiagram          zxDiag_reversed = zxDiag.reverse();
        qc::QuantumComputation qc_extracted_2  = qc::QuantumComputation(zxDiag_reversed.getNQubits());

        //Measurement measurement = Measurement(true); // FIXME:

        std::map<size_t, int> claimed_vertices;

        ExtractorParallel extractor1(qc_extracted, zxDiag, 0, &claimed_vertices, measurement);
        ExtractorParallel extractor2(qc_extracted_2, zxDiag_reversed, 1, &claimed_vertices, measurement);

        extractor1.other_extractor = &extractor2; // IMPROVE: needed?
        extractor2.other_extractor = &extractor1;

        std::cout << "Starting parallel extraction" << std::endl;
        #pragma omp parallel num_threads(2) shared(claimed_vertices)
        {
            if (omp_get_thread_num() == 0) {
                extractor1.extract();
                std::cout << "Extractor 1 finished" << std::endl;
            } else {
                extractor2.extract();
                std::cout << "Extractor 2 finished" << std::endl;
                //std::cout << qc_extracted_2 << std::endl;
                //qc_extracted_2.dump("H:/Uni/Masterarbeit/pyzx/thesis/extracted2.qasm");
            }
        }

        // TODO: Extract rest
        std::cout << "Finalizing extraction" << std::endl;
        extractor1.finalizeExtraction(extractor2.frontier);
        qc_extracted.dump("H:/Uni/Masterarbeit/pyzx/thesis/extracted2.qasm");

        // TODO: Combine diagrams
        std::cout << "Combining diagrams" << std::endl;
        qc_extracted_2.combine(qc_extracted);

        end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract", begin, end);

        std::cout << "Finished Circuit" << std::endl;
        std::cout << qc_extracted_2 << std::endl;
        //std::cout << "Circuit to extract:" << std::endl;
        //std::cout << qc << std::endl;
       
        qc_extracted_2.dump("H:/Uni/Masterarbeit/pyzx/thesis/extracted.qasm");
        std::cout << "Circuit " << circuitName << ":" << std::endl;

        measurement.printMeasurements(measurementGroup, circuitName);
    }
} // namespace zx