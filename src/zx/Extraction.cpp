#include "zx/Extraction.hpp"


#include "Definitions.hpp"
#include "Rational.hpp"
#include "Utils.hpp"
#include "ZXDiagram.hpp"
#include "Simplify.hpp"
#include "dd/Control.hpp"
#include "zx/FunctionalityConstruction.hpp"
#include "QuantumComputation.hpp"

//#include "EquivalenceCheckingManager.hpp"
//#include "gtest/gtest.h"

#include <algorithm>
#include <cstddef>
#include <optional>

#include <chrono>
#include <tuple>

/* struct VertexData {
        Col          col;
        Qubit        qubit;
        PiExpression phase;
        VertexType   type;
    }; */

using namespace dd;

namespace zx {

    void printVector(std::vector<size_t> vec) {
        for(auto const& v : vec) {
            std::cout << v << " ,";
        }
        std::cout << std::endl;
    }

    void printVector(std::vector<int> vec) {
        for(auto const& v : vec) {
            std::cout << v << " ,";
        }
        std::cout << std::endl;
    }

    void printVector(std::map<int, int> vec) {
        for(auto const& v : vec) {
            std::cout << v.first << " : " << v.second << " ,";
        }
        std::cout << std::endl;
    }

    void printVector(std::vector<bool> vec) {
        for(auto const& v : vec) {
            std::cout << v << " ,";
        }
        std::cout << std::endl;
    }

    void printVector(std::vector<std::pair<int, int>> vec) {
        for(auto const& v : vec) {
            std::cout << "(" << v.first << "," << v.second << ") ,";
        }
        std::cout << std::endl;
    }

    void printVector(std::set<std::pair<zx::Vertex, zx::Vertex>> vec) {
        for(auto const& v : vec) {
            std::cout << "(" << v.first << "," << v.second << ") ,";
        }
        std::cout << std::endl;
    }

    void printMatrix(gf2Mat matrix) {
        for(auto const& row : matrix) {
            printVector(row);
        }
        std::cout << std::endl;
    }
    




    /* initFrontier(zx::ZXDiagram& diag) {


    } */





    

    void extract(qc::QuantumComputation& circuit, zx::ZXDiagram& diag) {
        diag.toGraphlike(); // TODO: Not necessarily necessary

        //  TODO: Create a copy of diag

        // Create empty circuit
        
        
        // Initialise frontier
        auto inputs = diag.getInputs();
        auto outputs = diag.getOutputs();
        auto frontier = diag.getConnectedSet(outputs);

        //begin = std::chrono::steady_clock::now();
        //end = std::chrono::steady_clock::now();
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;

        
        std::cout << "Inputs:" << std::endl;
        printVector(inputs);

        std::cout << "Frontier:" << std::endl;
        printVector(frontier);

        std::cout << "Outputs:" << std::endl;
        printVector(outputs);
   
        
        processFrontier(diag, frontier, outputs, circuit);
        std::cout << "Current Circuit" << std::endl;
        std::cout << circuit << std::endl;

        int i = 1;

        while(updateFrontier(diag, frontier, outputs, circuit)) {
            std::cout << "Iteration " << i << std::endl;
            std::cout << "Current Circuit" << std::endl;
            std::cout << circuit << std::endl;
            i++;
        };

       

        /* for(auto const& v : frontier) {
            std::cout << "ASDASD" << std::endl;
            // If v connected to input by Hadamard
            auto incident = diag.incidentEdges(v);
            for(auto const &edge : incident) {
                zx::Vertex w = edge.to;
                if(diag.type(w) == zx::VertexType::Boundary) {
                    if(edge.type == zx::EdgeType::Hadamard) {
                        circuit.h(diag.qubit(v));
                        diag.removeEdge(v, w); 
                    }
                }
            } 
        } */

        // Reverse circuit
        

        // TODO: SWAP
        inputs = diag.getInputs();
        outputs = diag.getOutputs();
        //frontier = diag.getConnectedSet(outputs);

        std::cout << "Inputs:" << std::endl;
        printVector(inputs);

        std::cout << "Frontier:" << std::endl;
        printVector(frontier);

        std::cout << "Outputs:" << std::endl;
        printVector(outputs);

        std::cout << "EDGES: " << std::endl;
        for(  auto v : diag.getEdges() ) {
            std::cout << "(" << v.first << ", " << v.second << "), ";
        }
        std::cout << std::endl;


        //std::vector<std::pair<int, int>> swaps;
        std::map<int, int> swaps;
        bool leftover_swaps = false;

        for(size_t i = 0; i < frontier.size(); ++i) {
            zx::Vertex v = frontier[i];
            
            auto incident = diag.incidentEdges(v);


            for(auto &edge : incident) {
                zx::Vertex w = edge.to;
                auto it = std::find(inputs.begin(), inputs.end(), w);
                if( it != inputs.end()) {
                    // Check if edge to input is hadamard
                    if(edge.type == zx::EdgeType::Hadamard) {
                        circuit.h(diag.qubit(v));
                        //diag.removeEdge(v, w); 
                        edge.type = EdgeType::Simple;
                    }
                    size_t j = it - inputs.begin();
                    if( j != i) {
                        std::cout << "Found swap at " << i << " , " << j << std::endl;
                        leftover_swaps = true;
                    }
                    swaps[i] = j;
                    break;
                }
            } 
            //if(inputs[])
        }

        /* if(leftover_swaps) { // FIXME: There appears to be a bug...
            std::cout << "Creating swaps... " << std::endl;
            // Check for swaps
            for(auto s : permutation_as_swaps(swaps)) {
                circuit.swap(s.first, s.second);
            }
        } */
        

        circuit.reverse();

        return;
    }

    std::vector<std::pair<int, int>> permutation_as_swaps(std::map<int, int> perm) {
        std::cout << "Perm:" << std::endl;
        printVector(perm);

        std::vector<std::pair<int, int>> swaps;
        std::cout << "Size: " << perm.size() << std::endl;
        std::cout << "Size: " << perm.size() << std::endl;

        std::vector<int> l;
        for(size_t i = 0; i < perm.size(); ++i) {
            l.emplace_back(perm[i]);
        }

        std::cout << "l:" << std::endl;
        printVector(l);

        std::map<int, int> pinv; 
        for(auto i : perm) {
            pinv[i.second] = i.first;
        }

        std::cout << "pinv:" << std::endl;
        printVector(pinv);

        std::vector<int> linv;
        for(size_t i = 0; i < pinv.size(); ++i) {
            linv.emplace_back(pinv[i]);
        }

        std::cout << "linv:" << std::endl;
        printVector(linv);

        for(size_t i = 0; i < perm.size(); ++i) {
            if((size_t) l[i] == i) continue;
            int t1 = l[i];
            int t2 = linv[i];
            std::cout << "Adding swap gate at " << i << " , " << t2 << std::endl;
            swaps.emplace_back(std::pair<int,int>(i,t2));
            l[t2] = t1;
            linv[t1] = t2;
        }
        return swaps;
    }

    /* std::vector<std::pair<int, int>> gaussElimination(zx::gf2Mat& matrix) { // FIXME: Does not reduce fully
        int n = matrix.size();
        int m = matrix[0].size();
        std::cout << n << " x " << m << std::endl;
        std::vector<std::pair<int, int>> rowOperations;
        
        // For every row
        for(int i = 0; i < n; ++i){
            //std::cout << "Row: " << i << " ,Entry:" << matrix[i][i] << std::endl;
            // Entry is 0: Pivoting
            if(matrix[i][i] == false) {

                if(i < m-1) { // Special case: non-quadratic
                    for(int l = i; l < m; ++l) {
                        if(matrix[i][l] != false) {
                            // For every row
                            for(int j = 0; j < n; ++j){
                                if(j != i && matrix[j][l] != false) {
                                    // For every column
                                    for(int k = 0; k < m; ++k) {
                                        matrix[j][k] = matrix[j][k] ^ matrix[i][k];
                                    }
                                    rowOperations.emplace_back(std::pair<int, int>{i, j});
                                }
                            }
                        }
                    }
                    break;;
                }
                std::cout << "Trying to pivot!" << std::endl;
                for(int j = i+1; j < n; ++j) {
                    if(matrix[j][i] != false) {
                        // Swap row i with row j
                        for(int k = 0; k < m; ++k) {
                            std::cout << "Swapping rows " << i << " and " << j << std::endl;
                            int temp = matrix[i][k];
                            matrix[i][k] = matrix[j][k];
                            matrix[j][k] = temp;
                            // TODO: What to do here?
                        }
                        break;
                    }
                    else if(j == n-1) {
                        std::cout << "Giving up!" << std::endl;
                        return rowOperations;
                    }
                }
            }

            // For every row
            for(int j = 0; j < n; ++j){
                if(j != i && matrix[j][i] != false) {
                    // For every column
                    for(int k = 0; k < m; ++k) {
                        matrix[j][k] = matrix[j][k] ^ matrix[i][k];
                    }
                    rowOperations.emplace_back(std::pair<int, int>{i, j});
                }
            }
        }
        return rowOperations;
    } */

    std::vector<std::pair<int, int>> gaussElimination(zx::gf2Mat& matrix) {
        int n = matrix.size();
        int m = matrix[0].size();
        std::cout << n << " x " << m << std::endl;
        std::vector<std::pair<int, int>> rowOperations;
        
        // For every row
        for(int i = 0; i < n; ++i){
            //std::cout << "Row: " << i << " ,Entry:" << matrix[i][i] << std::endl;
            // Entry is 0: Pivoting
            if(matrix[i][i] == false) {

                if(i >= n) { // Special case: non-quadratic, then move further to the right to eliminate more entries
                    std::cout << "Non-quadratic matrix, moving to the right!" << std::endl;
                    for(int l = i; l < m; ++l) { // For every column to the right
                        if(matrix[i][l] != false) {
                            for(int j = 0; j < n; ++j){ // For every row
                                if(j != i && matrix[j][l] != false) {
                                    for(int k = 0; k < m; ++k) { // For every column
                                        matrix[j][k] = matrix[j][k] ^ matrix[i][k];
                                    }
                                    rowOperations.emplace_back(std::pair<int, int>{i, j});
                                }
                            }
                        }
                    }
                    break; // Breaks out of main loop; we reduced everything we could
                }

                //std::cout << "Current row is zero. Trying to find a non-zero row!" << std::endl;
                for(int j = i+1; j < n; ++j) { // For every row below the current
                    if(matrix[j][i] != false) {
                        // Add the other row onto this one
                        for(int k = 0; k < m; ++k) {
                            //std::cout << "Adding rows " << i << " and " << j << std::endl;
                            matrix[i][k] = matrix[j][k] ^ matrix[i][k];
                        }
                        rowOperations.emplace_back(std::pair<int, int>{j, i});
                        break;
                    }
                    else if(j == n-1) { // FIXME: Just continue...
                        std::cout << "Giving up!" << std::endl;
                        return rowOperations;
                    }
                }
            }

            // For every row
            for(int j = 0; j < n; ++j){
                if(j != i && matrix[j][i] != false) { // Add the current row onto it if it has a non-zero entry
                    // For every column
                    for(int k = 0; k < m; ++k) {
                        matrix[j][k] = matrix[j][k] ^ matrix[i][k];
                    }
                    rowOperations.emplace_back(std::pair<int, int>{i, j});
                }
            }
        }
        // TODO: Remove duplicates
        // Do (2,1) and (1,2) remove each other?
        return rowOperations;
    }



    // TODO: Move to ZXDiagram.hpp
    gf2Mat getAdjacencyMatrix(zx::ZXDiagram& diag, const std::vector<Vertex>& vertices_from, const std::vector<Vertex>& vertices_to) {
        gf2Mat adjMat{vertices_from.size(), gf2Vec(vertices_to.size(), false)};
        for(size_t i = 0; i < vertices_from.size(); ++i){
            for(size_t j = 0; j < vertices_to.size(); ++j){
                if(diag.connected(vertices_from[i], vertices_to[j])) {
                    adjMat[i][j] = true;
                }
                else {
                    adjMat[i][j] = false;
                }
            }
        }
        return adjMat;
    }

    // TODO: Move to ZXDiagram.hpp ?
    std::set<std::pair<Vertex, Vertex>> getEdgesBetweenVertices(zx::ZXDiagram& diag, std::vector<zx::Vertex> vertices) {
        std::set<std::pair<Vertex, Vertex>> ret;
        for(zx::Vertex v : vertices) {
            for(zx::Edge e : diag.incidentEdges(v)) {
                if( std::find(vertices.begin(), vertices.end(), e.to) != vertices.end() ) {
                    if(v < e.to) {
                        ret.insert( std::pair<Vertex, Vertex>{v, e.to} );
                    }
                    else {
                        ret.insert( std::pair<Vertex, Vertex>{e.to, v} );
                    }
                }
                
            }
        }
        return ret;
    }
    /* Alternative implementation. Seems to be about as fast...
    std::set<std::pair<Vertex, Vertex>> fe;
    for (auto v = frontier.begin(); v != frontier.end(); ++v) {
        for (auto w = v; ++w != frontier.end();) {
            if(diag.connected(*v, *w)) {
                fe.insert( std::pair<Vertex, Vertex>{*v, *w} );
            }
        }
    }
 */
    void processFrontier(zx::ZXDiagram& diag, std::vector<zx::Vertex>& frontier, std::vector<zx::Vertex>& outputs, qc::QuantumComputation& circuit) {
        std::cout << "VERTS: " << std::endl;
        for(  auto v : diag.getVertices() ) {
            std::cout << v.first << ", ";
        }
        std::cout << std::endl;

        std::cout << "EDGES: " << std::endl;
        for(  auto v : diag.getEdges() ) {
            std::cout << "(" << v.first << ", " << v.second << "), ";
        }
        std::cout << std::endl;

        for(auto const& v : frontier) {

            

            // If v connected to output by Hadamard
            // TODO: Improve to detect chains
            auto incident = diag.incidentEdges(v);
            for(auto edge : incident) {
                zx::Vertex w = edge.to;

                if( std::find(outputs.begin(), outputs.end(), w) == outputs.end() ) {
                    continue;
                }
                    
                if(edge.type == zx::EdgeType::Hadamard) {
                    circuit.h(diag.qubit(v));
                    edge.type = EdgeType::Simple;
                    //diag.removeEdge(v, w); // Remove edge, or remove Hadamard?
                }
                
            }

            // Add phase-gate at v with phase p
            if(!diag.phase(v).isZero()) {
                circuit.phase(diag.qubit(v), diag.phase(v).getConst().toDouble()); 
                diag.setPhase(v, PiExpression());
            }
            
        }
        // For every edge in frontier
        std::set<std::pair<Vertex, Vertex>> frontier_edges = getEdgesBetweenVertices(diag, frontier);
        std::cout << "Frontier Edges:" << std::endl;
        printVector(frontier_edges);

        for(std::pair<Vertex, Vertex> e : frontier_edges) {
            // Add cz-gate between v and w
            dd::Qubit qv = diag.qubit(e.first);
            circuit.z(diag.qubit(e.second), dd::Control{qv});

            // Remove edge between v and w
            diag.removeEdge(e.first, e.second);
        }        
    }

    bool updateFrontier(zx::ZXDiagram& diag, std::vector<zx::Vertex>& frontier, std::vector<zx::Vertex>& outputs, qc::QuantumComputation& circuit) {
        if(frontier.size() <= 0) return false;

        bool vertsRemaining = false;
        auto verts = diag.getVertices();

        std::vector<zx::Vertex> temp;

        // Check for remaining vertices
        for(auto v : verts) {
            if(diag.type(v.first) != zx::VertexType::Boundary && std::find(frontier.begin(), frontier.end(), v.first) == frontier.end()) {
                vertsRemaining = true;
                temp.emplace_back(v.first);
            }
        }

        if(!vertsRemaining) {
            std::cout << "No verts remaining! " << std::endl;
            return false;
        }

        std::cout << "Frontier:" << std::endl;
        printVector(frontier);

        std::cout << "Verts remaining:" << std::endl;
        printVector(temp);

        // Get neighbours of frontier
        auto neighbours = diag.getConnectedSet(frontier, outputs);

        std::cout << "Neighbours:" << std::endl;
        printVector(neighbours);

        // Get biadjacency matrix of frontier/neighbours
        auto adjMatrix = getAdjacencyMatrix(diag, frontier, neighbours);

        std::cout << "Adjacency Matrix:" << std::endl;
        printMatrix(adjMatrix);

        // Gauss reduction on biadjacency matrix
        std::vector<std::pair<int, int>> rowOperations = gaussElimination(adjMatrix);
        std::cout << "After Gauss Elimination:" << std::endl;
        printMatrix(adjMatrix);
        std::cout << "Row Operations:" << std::endl;
        printVector(rowOperations);

        std::vector<zx::Vertex> ws;
        std::vector<zx::Vertex> ws_neighbours;
        for(size_t i = 0; i < adjMatrix.size(); ++i) {
            int sum = 0, nonZero = 0;
            for(size_t j = 0; j < adjMatrix[i].size(); ++j) {
                sum += adjMatrix[i][j];
                if(adjMatrix[i][j]) {
                    nonZero = j;
                }
            }
            if(sum == 1) {
                // Set w to vertex corresponding to nonZero
                int w = neighbours[nonZero];
                ws.emplace_back(w); // FIXME: Check for duplicates (?)
                ws_neighbours.emplace_back(frontier[i]); // TODO: Make this into a map frontier[i] --> w
            }
        }
        std::cout << "Vector ws:" << std::endl;
        printVector(ws);

        std::cout << "Vector ws_neighbours:" << std::endl;
        printVector(ws_neighbours);

        if(ws.size() <= 0) {
            std::cout << "Ws is 0" << std::endl;
            // TODO:
            return false;
            while(true) {
                auto neighbours = diag.getConnectedSet(frontier, outputs);
                if(neighbours.size() <= 0) break;


            }
            

            processFrontier(diag, frontier, outputs, circuit);
            return true; // FIXME: ?
        }

        //adjMatrix = getAdjacencyMatrix(diag, frontier, ws);
        //std::cout << "New Adjacency Matrix:" << std::endl;
        //printMatrix(adjMatrix);

        // Extract CNOTs
        for(auto r : rowOperations) {
            dd::Qubit control = diag.qubit(r.second); // From row operation: second = second + first
            circuit.x(diag.qubit(r.first), dd::Control{control}); // FIXME: Which order is correct?

            // TODO: Update diag based on row operation // FIXME: Which order is correct?
            for(auto v : diag.getNeighbourVertices(r.first)) {

                if( std::find(outputs.begin(), outputs.end(), v) != outputs.end() ) {
                    continue;
                }

                if( std::find(frontier.begin(), frontier.end(), v) != frontier.end() ) {
                    continue;
                }

                if(diag.connected(v, r.second)) {
                    diag.removeEdge(v, r.second);
                }
                else {
                    if(diag.type(v) == zx::VertexType::Boundary) { // v is an input
                        auto new_v = diag.insertIdentity(r.first, v);
                        if(new_v) {
                            diag.addEdge(r.second, *new_v, zx::EdgeType::Hadamard);
                        }
                    }
                    else {
                        diag.addEdge(r.second, v, zx::EdgeType::Hadamard);
                    }
                }
            }
        }
        
        for(size_t i = 0; i < ws.size(); ++i) {
            zx::Vertex v = ws_neighbours[i];
            zx::Vertex w = ws[i];
            circuit.h(diag.qubit(v));
            circuit.phase(diag.qubit(v), diag.phase(w).getConst().toDouble()); 
            // TODO: 
            frontier.erase(std::remove(frontier.begin(), frontier.end(), v), frontier.end());  // TODO: Should really have been a set to begin with
            frontier.emplace_back(w);
            diag.removeVertex(v);
            
        }

        /* std::cout << "Frontier:" << std::endl;
        printVector(frontier);

        std::set<std::pair<Vertex, Vertex>> ws_edges = getEdgesBetweenVertices(diag, ws);
        std::cout << "ws_edges:" << std::endl;
        printVector(ws_edges);

        for(std::pair<Vertex, Vertex> e : ws_edges) {
            dd::Qubit qv = diag.qubit(e.second);
            circuit.z(diag.qubit(e.first), dd::Control{qv});
            diag.removeEdge(e.first, e.second);
        } */
        processFrontier(diag, frontier, outputs, circuit);
        return true;
    }

    void testExtraction() {
        std::cout << "Setting up...\n";
        qc::QuantumComputation qc{};
        qc.addQubitRegister(4U);
        /* qc.z(0, 3_pc);
        qc.z(1, 2_pc);
        qc.h(3);
        qc.phase(1, dd::PI_2);
        qc.phase(0, dd::PI_2);
        
        qc.x(0); */
        qc.x(1);
        qc.h(3); 
        qc.x(3, 2_pc); 
        qc.t(0);
        qc.t(1);
        qc.t(2);
        /*qc.tdag(3);
         qc.x(1, 0_pc);
        qc.x(3, 2_pc);
        qc.x(0, 3_pc);
        qc.x(2, 1_pc);
        qc.x(1, 0_pc);
        qc.x(3, 2_pc);
        qc.tdag(0);
        qc.tdag(1);
        qc.tdag(2);
        qc.t(3);
        qc.x(1, 0_pc);
        qc.x(3, 2_pc);
        qc.s(3);
        qc.x(0, 3_pc);
        qc.h(3); */
        qc.dump("D:/Uni/Masterarbeit/qcec/original.qasm");

        std::cout << "Circuit to extract:" << std::endl;
        std::cout << qc << std::endl;
        zx::ZXDiagram zxDiag = zx::FunctionalityConstruction::buildFunctionality(&qc);
        //zx::fullReduce(zxDiag);
        zx::cliffordSimp(zxDiag);
        qc::QuantumComputation qc_extracted = qc::QuantumComputation(zxDiag.getNQubits());
        extract(qc_extracted, zxDiag);
        std::cout << "Finished Circuit" << std::endl;
        std::cout << qc_extracted << std::endl;
        qc_extracted.dump("D:/Uni/Masterarbeit/qcec/extracted.qasm");

        /* ec::Configuration config{};
        config.functionality.traceThreshold    = 1e-2;
        config.execution.runAlternatingChecker = true;
        ec::EquivalenceCheckingManager ecm(qc, qc_extracted, config);
        ecm.run();
        std::cout << ecm << std::endl;
        EXPECT_EQ(ecm.equivalence(), ec::EquivalenceCriterion::Equivalent); */

    }
}



    

   



