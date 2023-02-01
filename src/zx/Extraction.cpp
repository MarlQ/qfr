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
#include <map>
#include <fstream>
#define DEBUG false

/* struct VertexData {
        Col          col;
        Qubit        qubit;
        PiExpression phase;
        VertexType   type;
    }; */

using namespace dd;

namespace zx {

    /* template <typename T>
    void printVector(const std::vector<T>& vec) {
        for (auto i = vec.begin(); i != vec.end(); ++i)
            if(DEBUG)std::cout << *i << ' ';
        if(DEBUG)std::cout << std:: endl;
    } */

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


    bool contains(std::vector<size_t> vec, zx::Vertex val) {        // IMPROVE: isIn is already implemented in ZXDiagram
        return std::find(vec.begin(), vec.end(), val) != vec.end();
    }

    bool contains(std::map<zx::Qubit, zx::Vertex> vec, zx::Vertex val) {
        return std::find_if(vec.begin(), vec.end(), [val](const auto& p) { return p.second == val; }) != vec.end();
    }

    std::map<std::string, double> times = {
        {"extract:extractRZ_CZ:PhaseGates", 0.0},
        {"extract:extractRZ_CZ:CZGates", 0.0},
        {"extract:extractRZ_CZ", 0.0},
        {"extract:extractCNOT", 0.0},
        {"extract:processFrontier", 0.0},
        {"extract:extractCNOT:gaussElimination", 0.0},
        {"extract:extractCNOT:CNOTFromOperations", 0.0},
        {"buildFunctionality", 0.0},
        {"interiorCliffordSimp", 0.0},
        {"extract", 0.0},
        {"extract:extractCNOT:YZSpiders", 0.0},
        {"extract:extractCNOT:biadjacencyMatrix", 0.0},
        {"extract:extractCNOT:check", 0.0}
    };

    void addMeasurement(std::string name, std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end) {
        double duration = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count();
        times[name] += duration;
        //std::cout << name << " = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
    }

    void printMeasurements(std::string iteration, std::string filename) {
        for(auto t : times) {
            std::cout << t.first << " = " << t.second / 1000000.0 << "[ms]" << std::endl;
        }

        std::ofstream file("H:/Uni/Masterarbeit/measurements.csv", std::ios::app);
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::stringstream date;
        date << std::put_time(std::localtime(&now_c), "%Y-%m-%d %X");
        
        file << iteration << "," << filename << "," << date.str();
        for(const auto& [key, value] : times) {
            file << "," << value / 1000000.0;
        }
        file << std::endl;
        file.close();
    }

    
    // IMPROVE: Move to ZXDiagram.hpp ?
    std::set<std::pair<Vertex, Vertex>> getEdgesBetweenVertices(zx::ZXDiagram& diag, std::vector<zx::Vertex> vertices) {
        std::set<std::pair<Vertex, Vertex>> ret;
        for(zx::Vertex v : vertices) {
            for(zx::Edge e : diag.incidentEdges(v)) {
                if( contains(vertices, e.to) ) {
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

    std::set<std::pair<Vertex, Vertex>> getEdgesBetweenVertices(zx::ZXDiagram& diag, std::map<zx::Qubit, zx::Vertex> vertices) {
        std::set<std::pair<Vertex, Vertex>> ret;
        for(auto v : vertices) {
            for(zx::Edge e : diag.incidentEdges(v.second)) {
                if( contains(vertices, e.to) ) {
                    if(v.second < e.to) {
                        ret.insert( std::pair<Vertex, Vertex>{v.second, e.to} );
                    }
                    else {
                        ret.insert( std::pair<Vertex, Vertex>{e.to, v.second} );
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

    std::vector<Vertex> mapToVector(std::map<zx::Qubit, zx::Vertex>& vertices) {
        std::vector<Vertex> retVerts;

        for( std::map<zx::Qubit, zx::Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it ) {
            retVerts.push_back( it->second );
         }
        return retVerts;
    }


    void extractRZ_CZ(zx::ZXDiagram& diag, std::map<zx::Qubit, zx::Vertex>& frontier, qc::QuantumComputation& circuit) {
        if(DEBUG)std::cout << "Extracting RZ and CZ gates..." << std::endl;

        auto begin = std::chrono::steady_clock::now();
        for(auto v : frontier) {
            // Extract RZ: Add phase-gate at v with phase p
            if(!diag.phase(v.second).isZero()) {
                if(DEBUG)std::cout << "Adding phase gate at " << v.first << " with phase " << diag.phase(v.second).getConst().toDouble() << std::endl;
                if(!diag.phase(v.second).isConstant()) {
                    if(DEBUG)std::cout << "Error: phase is not constant!" << std::endl;
                    exit(0);
                }
                circuit.phase(v.first, diag.phase(v.second).getConst().toDouble()); 
                diag.setPhase(v.second, PiExpression());
            }
        }
        auto end = std::chrono::steady_clock::now();
        addMeasurement("extract:extractRZ_CZ:PhaseGates", begin, end);
        
        begin = std::chrono::steady_clock::now();
        // Extract connected frontier spiders as CZ gates
        //std::set<std::pair<Vertex, Vertex>> frontier_edges = getEdgesBetweenVertices(diag, frontier);
        //if(DEBUG)std::cout << "Frontier Edges:" << std::endl;
        //if(DEBUG)printVector(frontier_edges);
        
        for(auto v : frontier) { // IMPROVE: Can this be within the same for loop as the prior?
            for(zx::Edge e : diag.incidentEdges(v.second)) {
                auto w = e.to;
                auto it = std::find_if(frontier.begin(), frontier.end(), [w](const auto& p) { return p.second == w; });
                if( it != frontier.end() ) {
                    dd::Qubit qw = it->first;
                    if(DEBUG)std::cout << "Adding CZ gate at " << v.first << "/" << it->first << std::endl;
                   
                    circuit.z(v.first, dd::Control{qw});

                    // Remove edge between v and w
                    diag.removeEdge(v.second, w);
                }
            }
        }
        end = std::chrono::steady_clock::now();
        addMeasurement("extract:extractRZ_CZ:CZGates", begin, end);
    }

    std::map<zx::Qubit, zx::Vertex> initFrontier(zx::ZXDiagram& diag) {
        std::map<zx::Qubit, zx::Vertex> frontier;
        auto outputs = diag.getOutputs();
        for(size_t i = 0; i < outputs.size(); ++i) {
            auto v = diag.getNeighbourVertices(outputs[i])[0];
            if( !diag.isInput(v) ) {
                frontier[i] = v;
            }
        }

        return frontier;
    }

    void extract(qc::QuantumComputation& circuit, zx::ZXDiagram& diag) {
        diag.toGraphlike(); // IMPROVE: Not necessarily necessary

        //  IMPROVE: Create a copy of diag

        // Create empty circuit
        
        
        // Initialise frontier
        auto inputs = diag.getInputs();
        auto outputs = diag.getOutputs();
        std::map<zx::Qubit, zx::Vertex> frontier = initFrontier(diag);

        

        
        if(DEBUG)std::cout << "Inputs:" << std::endl;
        if(DEBUG)printVector(inputs);

        if(DEBUG)std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        if(DEBUG)std::cout << "Outputs:" << std::endl;
        if(DEBUG)printVector(outputs);

        for(auto v : frontier) { // IMPROVE: Iterate over outputs
            auto incident = diag.incidentEdges(v.second);
            for(auto edge : incident) {
                zx::Vertex w = edge.to;

                if( !diag.isOutput(w) ) {
                    continue;
                }
                    
                if(edge.type == zx::EdgeType::Hadamard) {
                    if(DEBUG)std::cout << "Adding Hadamard at " << v.first << std::endl;
                    circuit.h(v.first);
                    //edge.type = zx::EdgeType::Simple;
                    diag.setEdgeType(v.second, w, zx::EdgeType::Simple);
                    //edge.toggle();
                    // toggle corresponding edge in other direction
                    //diag.getEdgePtr(edge.to, v)->toggle();
                }
                
            }
        } 
        
        if(DEBUG)std::cout << "Current Circuit" << std::endl;
        if(DEBUG)std::cout << circuit << std::endl;

        int i = 1;

        while(frontier.size() > 0) {
            auto begin = std::chrono::steady_clock::now();
            extractRZ_CZ(diag, frontier, circuit);
            auto end = std::chrono::steady_clock::now();
            addMeasurement("extract:extractRZ_CZ", begin, end);

            begin = std::chrono::steady_clock::now();
            extractCNOT(diag, frontier, outputs, circuit);
            end = std::chrono::steady_clock::now();
            addMeasurement("extract:extractCNOT", begin, end);

            begin = std::chrono::steady_clock::now();
            processFrontier(diag, frontier, outputs, circuit);
            end = std::chrono::steady_clock::now();
            addMeasurement("extract:processFrontier", begin, end);
            if(DEBUG)std::cout << "Iteration " << i << std::endl;
            //std::cout << "Current Circuit" << std::endl;
            //std::cout << circuit << std::endl;
            i++;
            //if(i == 6) exit(0);
        };

        // Reverse circuit
        if(DEBUG)std::cout << "Finished extraction. Reversing circuit and finding swaps..." << std::endl;

        inputs = diag.getInputs();
        outputs = diag.getOutputs();
        //frontier = diag.getConnectedSet(outputs);

        if(DEBUG)std::cout << "Inputs:" << std::endl;
        if(DEBUG)printVector(inputs);

        if(DEBUG)std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        if(DEBUG)std::cout << "Outputs:" << std::endl;
        if(DEBUG)printVector(outputs);

        if(DEBUG)std::cout << "EDGES: " << std::endl;
        for(  auto v : diag.getEdges() ) {
            if(DEBUG)std::cout << "(" << v.first << ", " << v.second << "), ";
        }
        if(DEBUG)std::cout << std::endl;

        std::map<int, int> swaps;
        bool leftover_swaps = false;
        

        for(int q = 0; q < outputs.size(); ++q) {
            
            auto incident = diag.incidentEdges(outputs[q]);

            for(auto &edge : incident) {
                zx::Vertex w = edge.to;
                auto it = std::find(inputs.begin(), inputs.end(), w);
                if( it != inputs.end()) {
                    // Check if edge to input is hadamard
                    if(edge.type == zx::EdgeType::Hadamard) {
                        circuit.h(q);
                        //diag.removeEdge(v, w); 
                        diag.setEdgeType(outputs[q], edge.to, EdgeType::Simple);
                        //edge.type = EdgeType::Simple;
                    }
                    //size_t j = it - inputs.begin();
                    auto qw = (zx::Qubit) (it - inputs.begin()); // FIXME: Not sure if this is correct
                    if( qw != q) {
                        if(DEBUG)std::cout << "Found swap at " << q << " , " << qw << std::endl;
                        leftover_swaps = true;
                    }
                    swaps[ q ] = qw;
                    break;
                }
            } 
        }

        if(leftover_swaps) {
            if(DEBUG)std::cout << "Creating swaps... " << std::endl;
            // Check for swaps
            for(auto s : permutation_as_swaps(swaps)) {
                circuit.swap(s.first, s.second);
            }
        }
        

        circuit.reverse();

        return;
    }

    std::vector<std::pair<int, int>> permutation_as_swaps(std::map<int, int> perm) {
        if(DEBUG)std::cout << "Perm:" << std::endl;
        if(DEBUG)printVector(perm);

        std::vector<std::pair<int, int>> swaps;
        if(DEBUG)std::cout << "Size: " << perm.size() << std::endl;
        if(DEBUG)std::cout << "Size: " << perm.size() << std::endl;

        std::vector<int> l;
        for(size_t i = 0; i < perm.size(); ++i) {
            l.emplace_back(perm[i]);
        }

        if(DEBUG)std::cout << "l:" << std::endl;
        if(DEBUG)printVector(l);

        std::map<int, int> pinv; 
        for(auto i : perm) {
            pinv[i.second] = i.first;
        }

        if(DEBUG)std::cout << "pinv:" << std::endl;
        if(DEBUG)printVector(pinv);

        std::vector<int> linv;
        for(size_t i = 0; i < pinv.size(); ++i) {
            linv.emplace_back(pinv[i]);
        }

        if(DEBUG)std::cout << "linv:" << std::endl;
        if(DEBUG)printVector(linv);

        for(size_t i = 0; i < perm.size(); ++i) {
            if((size_t) l[i] == i) continue;
            int t1 = l[i];
            int t2 = linv[i];
            if(DEBUG)std::cout << "Adding swap gate at " << i << " , " << t2 << std::endl;
            swaps.emplace_back(std::pair<int,int>(i,t2));
            l[t2] = t1;
            linv[t1] = t2;
        }
        return swaps;
    }

     std::vector<std::pair<zx::Qubit, zx::Qubit>> gaussEliminationOther(zx::gf2Mat& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();
        int lead = 0;
        if(DEBUG)std::cout << rows << " x " << cols << std::endl;
        std::vector<std::pair<zx::Qubit, zx::Qubit>> rowOperations;
        for (int r = 0; r < rows; r++) {
            if (cols <= lead)
                return rowOperations;
            int i = r;
            while (matrix[i][lead] == 0) {
                i++;
                if (i == rows) {
                    i = r;
                    lead++;
                    if (cols == lead)
                        return rowOperations;
                }
            }
            /* if (i != r) {
                for (int j = 0; j < cols; j++) {
                    matrix[i][j] = matrix[i][j] ^ matrix[r][j];
                    matrix[r][j] = matrix[r][j] ^ matrix[i][j];
                    matrix[i][j] = matrix[i][j] ^ matrix[r][j];          
                }
                rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{r, i});
                rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{i, r});
                rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{r, i});
            } */
            for (int j = 0; j < rows; j++) {
                if (j != i) { // FIXME: Replaced r with i
                    if (matrix[j][lead]) {
                            for (int k = 0; k < cols; k++) { // IMPROVE: k = lead ?
                                matrix[j][k] = matrix[j][k] ^ matrix[i][k]; // FIXME: Replaced r with i
                            }
                            // Check if the row accidentally got set to 0
                            bool allZero = true;
                            for (int k = 0; k < cols; k++) {
                                if(matrix[j][k]) {
                                    allZero = false;
                                    break;
                                }
                            }
                            // Undo row addition
                            if(allZero) {
                                for (int k = 0; k < cols; k++) { // IMPROVE: k = lead ?
                                    matrix[j][k] = matrix[j][k] ^ matrix[i][k]; // FIXME: Replaced r with i
                                }
                            }
                            else {
                                rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{i, j});
                            }
                        }
                }
            }
            lead++;
        }
        return rowOperations;
    }

    // FIXME: Currently generates way more gates than it has to (because of echelon form)
    std::vector<std::pair<zx::Qubit, zx::Qubit>> gaussElimination(zx::gf2Mat& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();
        if(DEBUG)std::cout << rows << " x " << cols << std::endl;
        std::vector<std::pair<zx::Qubit, zx::Qubit>> rowOperations;
        
        // For every row
        for(int i = 0; i < rows; ++i){
            //std::cout << "Row: " << i << " ,Entry:" << matrix[i][i] << std::endl;
            // Entry is 0: Pivoting
            if(!matrix[i][i]) {

                if(i >= rows) { // Special case: non-quadratic, then move further to the right to eliminate more entries
                    if(DEBUG)std::cout << "Non-quadratic matrix, moving to the right!" << std::endl;
                    for(int l = i; l < cols; ++l) { // For every column to the right
                        if(matrix[i][l] != false) {
                            for(int j = 0; j < rows; ++j){ // For every row
                                if(j != i && matrix[j][l] != false) {
                                    for(int k = 0; k < cols; ++k) { // For every column
                                        matrix[j][k] = matrix[j][k] ^ matrix[i][k];
                                    }
                                    rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{i, j}); // FIXME:
                                }
                            }
                        }
                    }
                    break; // Breaks out of main loop; we reduced everything we could
                }

                //std::cout << "Current row is zero. Trying to find a non-zero row!" << std::endl;
                for(int j = i+1; j < rows; ++j) { // For every row below the current
                    if(matrix[j][i]) {
                        // Add the other row onto this one
                        for(int k = 0; k < cols; ++k) {
                            //std::cout << "Adding rows " << i << " and " << j << std::endl;
                            matrix[i][k] = matrix[j][k] ^ matrix[i][k];
                        }
                        rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{j, i}); // FIXME:
                        break;
                    }
                    else if(j == rows-1) {
                        /* if(DEBUG)std::cout << "Giving up!" << std::endl;
                        return rowOperations; */
                        goto nextRow; // Not amazing, but faster than using functions
                    }
                }
            }

            // For every row
            for(int j = 0; j < rows; ++j){
                if(j != i && matrix[j][i]) { // Add the current row onto it if it has a non-zero entry
                    // For every column
                    for(int k = 0; k < cols; ++k) {
                        matrix[j][k] = matrix[j][k] ^ matrix[i][k];
                    }
                    rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{i, j}); // FIXME:
                }
            }
            nextRow: ;

        }
        // TODO: Remove duplicates
        // Do (2,1) and (1,2) remove each other?
        return rowOperations;
    }

    void row_add(zx::gf2Mat& matrix, int r0, int r1, std::vector<std::pair<zx::Qubit, zx::Qubit>> &rowOperations) { // CHECK: Should rowOperations be a set, instead of a vector to avoid duplicates?
        int cols = matrix[0].size();
        for(int k = 0; k < cols; ++k) {
            matrix[r1][k] = matrix[r1][k] ^ matrix[r0][k];
        }
        rowOperations.emplace_back(std::pair<zx::Qubit, zx::Qubit>{r0, r1});
    }

    std::vector<std::pair<zx::Qubit, zx::Qubit>> gaussEliminationAlt(zx::gf2Mat& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();
        int blocksize = 6;
        int pivot_row = 0;
        std::vector<std::pair<zx::Qubit, zx::Qubit>> rowOperations;
        std::vector<int> pivot_cols;

        for(int sec = 0; sec < std::ceil(cols / (double)blocksize); ++sec) {
                int i0 = sec * blocksize;
                int i1 = std::min(cols, (sec+1) * blocksize);

                std::unordered_map<std::vector<bool>, int> chunks;
                for(int r = pivot_row; r < rows; ++r) {
                    std::vector<bool> t;
                    t.reserve(i1 - i0);
                    std::copy(matrix[r].begin() + i0, matrix[r].begin() + i1, std::back_inserter(t));
                    if (std::none_of(t.begin(), t.end(), [](bool b) { return b;})) {
                        continue;
                    }
                    auto it = chunks.find(t);
                    if (it != chunks.end()) {
                        row_add(matrix, chunks[t], r, rowOperations);
                    }
                    else {
                        chunks[t] = r;
                    }
                }
                int p = i0;
                while(p < i1) {
                    for(int r0 = pivot_row; r0 < rows; ++r0) {
                        if(matrix[r0][p]) {
                            if(r0 != pivot_row) {
                                row_add(matrix, r0, pivot_row, rowOperations);
                            }

                            for(int r1 = pivot_row + 1; r1 < rows; ++r1) {
                                if(pivot_row != r1 && matrix[r1][p]) {
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
            int i1 = std::min(cols, (sec+1) * blocksize);

            std::unordered_map<std::vector<bool>, int> chunks;
            for(int r = pivot_row - 1; r >= 0; r--) {
                std::vector<bool> t;
                t.reserve(i1 - i0);
                std::copy(matrix[r].begin() + i0, matrix[r].begin() + i1, std::back_inserter(t));
                if (std::none_of(t.begin(), t.end(), [](bool b) { return b;})) {
                    continue;
                }
                auto it = chunks.find(t);
                if (it != chunks.end()) {
                    row_add(matrix, chunks[t], r, rowOperations);
                }
                else {
                    chunks[t] = r;
                }
            }

            while(pivot_cols1.size() != 0 && i0 <= pivot_cols1.back() < i1) {
                int pcol = pivot_cols1.back();
                pivot_cols1.pop_back();
                for(int r = 0; r < pivot_row; ++r) {
                    if(matrix[r][pcol]) {
                        row_add(matrix, pivot_row, r, rowOperations);
                    }
                }
                pivot_row--;
            }
        }
        return rowOperations;
    }


    // TODO: Move to ZXDiagram.hpp
    gf2Mat getAdjacencyMatrix(zx::ZXDiagram& diag, const std::map<zx::Qubit, zx::Vertex>& vertices_from, const std::vector<Vertex>& vertices_to) {
        if(DEBUG)std::cout << " Matrix: " << vertices_from.size() << " x " << vertices_to.size() << std::endl;
        gf2Mat adjMat{vertices_from.size(), gf2Vec(vertices_to.size(), false)};
        for(size_t i = 0; i < vertices_from.size(); ++i){
            for(size_t j = 0; j < vertices_to.size(); ++j){
                if(diag.connected(vertices_from.at(i), vertices_to[j])) {
                    adjMat[i][j] = true;
                }
                else {
                    adjMat[i][j] = false;
                }
            }
        }
        return adjMat;
    }

    gf2Mat getAdjacencyMatrix(zx::ZXDiagram& diag, const std::vector<zx::Vertex>& vertices_from, const std::vector<Vertex>& vertices_to) {
        if(DEBUG)std::cout << " Matrix: " << vertices_from.size() << " x " << vertices_to.size() << std::endl;
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

    

    void processFrontier(zx::ZXDiagram& diag, std::map<zx::Qubit, zx::Vertex>& frontier, std::vector<zx::Vertex>& outputs, qc::QuantumComputation& circuit) {
        if(DEBUG)std::cout << "Processing Frontier... " << std::endl;
        auto inputs = diag.getInputs();
        if(DEBUG)std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);
        /* if(DEBUG)std::cout << "Inputs:" << std::endl;
        printVector(inputs);

        if(DEBUG)std::cout << "Frontier:" << std::endl;
        printVector(frontier);

        if(DEBUG)std::cout << "Outputs:" << std::endl;
        printVector(outputs);

        if(DEBUG)std::cout << "VERTS: " << std::endl;
        for(  auto v : diag.getVertices() ) {
            if(DEBUG)std::cout << v.first << ", ";
        }
        if(DEBUG)std::cout << std::endl;

        if(DEBUG)std::cout << "EDGES: " << std::endl;
        for(  auto v : diag.getEdges() ) {
            if(DEBUG)std::cout << "(" << v.first << ", " << v.second << "), ";
        }
        if(DEBUG)std::cout << std::endl; */

        std::map<zx::Qubit, int> new_frontier;

        for(auto const& v : frontier) {
            if(DEBUG)std::cout << "Vertex: " << v.second << std::endl;
            auto current_neighbours = diag.getNeighbourVertices(v.second);
            if(DEBUG)std::cout << "Neighb.:" << std::endl;
            if(DEBUG)printVector(current_neighbours);
            if(current_neighbours.size() > 2 || contains(inputs, v.second)) {
                continue; // Process later
            }
            zx::Vertex output;

            // Eliminate Hadamard chains
            bool uneven_hadamard = false;
            zx::Vertex previous_vertex;
            zx::Vertex current_vertex = v.second;
            std::vector<zx::Vertex> chain;

            for(zx::Vertex n : current_neighbours) {
                if(diag.isOutput(n)) {
                    previous_vertex = n;
                    output = n;
                    //std::cout << "Output: " << output << std::endl;
                    break;
                }
            }

            while(true) {
                zx::Vertex next_vertex;
                //std::cout << "Current Vertex: " << current_vertex << std::endl;
                //std::cout << "Previous Vertex: " << previous_vertex << std::endl;
                
                for(auto n : current_neighbours) {
                    if(n != previous_vertex) {
                        next_vertex = n;
                    }
                }
                
                /* if(next_vertex == -1) {
                    if(DEBUG)std::cout << "ERROR: No next vertex!" << std::endl;
                    break;
                } */
                //std::cout << "Next Vertex: " << next_vertex<< std::endl;
                
                chain.emplace_back(next_vertex);

                auto edge = diag.getEdge(current_vertex, next_vertex);

                if(!edge) return;

                if(edge->type == EdgeType::Hadamard) {
                    uneven_hadamard = !uneven_hadamard;
                }
                if(!diag.phase(next_vertex).isZero()) { // IMPROVE: Should this really break?
                    break;
                }
                if(diag.getNeighbourVertices(next_vertex).size() > 2 || contains(inputs, next_vertex)) {
                    break;
                }

                previous_vertex = current_vertex;
                current_vertex = next_vertex;
                current_neighbours = diag.getNeighbourVertices(current_vertex);
            }

            if(uneven_hadamard) {
                circuit.h(v.first);
                if(DEBUG)std::cout << "Adding Hadamard at " << v.first << std::endl;
                //diag.setEdgeType(v.second, chain[0], zx::EdgeType::Simple); // CHECK: Is this a good way to do this?
            }
            if(DEBUG)std::cout << "Chain found:" << std::endl;
            if(DEBUG)printVector(chain);

            for(unsigned int i = 0; i < chain.size() - 1; i++) {
                if(DEBUG)std::cout << "Removing Vertex: " << chain[i] << std::endl;
                diag.removeVertex(chain[i]);
            }

            auto edgeType = diag.getEdge(v.second,output)->type; //CHECK: isn't this always Simple, as we removed hadamards at the start?
            auto last_in_chain = chain[chain.size() - 1];
            //auto v_qubit = v.first;
            if(DEBUG)std::cout << "Removing Vertex: " << v.second << std::endl;
            diag.removeVertex(v.second);
            if(DEBUG)std::cout << "Adding Edge: (" << last_in_chain << "," << output << ")" << std::endl;
            diag.addEdge(last_in_chain, output, edgeType);
            
            if(!contains(inputs, last_in_chain)) {
                new_frontier[ v.first ] = (int) last_in_chain; // FIXME: Potentially problematic cast
            }
            else {
                new_frontier[ v.first ] = -1;
            }
            //std::cout << "Frontier Changes:" << std::endl;
            if(DEBUG)printVector(new_frontier);

        }
        if(DEBUG)std::cout << "Old Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);
        
        if(new_frontier.size() > 0) {
            if(DEBUG)std::cout << "Frontier Changes:" << std::endl;
            if(DEBUG)printVector(new_frontier);
            for(auto entry : new_frontier) {
                //frontier.erase(std::remove(frontier.begin(), frontier.end(), entry.first), frontier.end());
                frontier.erase( entry.first );
                if(entry.second != -1) {
                    frontier[ entry.first ] = (zx::Vertex) entry.second;
                }
            }
        }     
        if(DEBUG)std::cout << "New Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);
        
    }

    std::vector<zx::Vertex> get_frontier_neighbors(zx::ZXDiagram& diag, std::map<zx::Qubit, zx::Vertex>& frontier) {
        std::vector<Vertex> frontier_neighbors;
        auto outputs = diag.getOutputs();

        for(auto v : frontier) {
            auto v_neighbours = diag.getNeighbourVertices(v.second);
            for(auto w: v_neighbours) {
                if( !contains(outputs, w) && !contains(frontier_neighbors, w) ) {
                    frontier_neighbors.emplace_back(w);
                }
            }
        }

        return frontier_neighbors;
    }

    void extractCNOT(zx::ZXDiagram& diag, std::map<zx::Qubit, zx::Vertex>& frontier, std::vector<zx::Vertex>& outputs, qc::QuantumComputation& circuit) {
        if(DEBUG)std::cout << "Is CNOT extraction necessary?" << std::endl;
        //bool vertsRemaining = false;
        auto verts = diag.getVertices();

        std::vector<zx::Vertex> temp;

        // Check for remaining vertices
        auto begin = std::chrono::steady_clock::now();
        for(auto v : frontier) { // FIXME: Is this correct?
            auto neighbours = diag.getNeighbourVertices(v.second);
            if(DEBUG)std::cout << "Neighbours of " << v.first << " :" << std::endl;
            if(DEBUG)printVector(neighbours);
            if(neighbours.size() <= 2) {
                if(DEBUG)std::cout << "No need for CNOT extraction. " << std::endl;
                return;
            }
        }
        auto end = std::chrono::steady_clock::now();
        addMeasurement("extract:extractCNOT:check", begin, end);
        if(DEBUG)std::cout << "Frontier CNOT extraction... " << std::endl;
/*         if(DEBUG)std::cout << "Current Circuit" << std::endl;
        if(DEBUG)std::cout << circuit << std::endl; */

        // OLD: Check for remaining vertices
        /* for(auto v : verts) {
            if(diag.type(v.first) != zx::VertexType::Boundary && std::find(frontier.begin(), frontier.end(), v.first) == frontier.end()) {
                vertsRemaining = true;
                temp.emplace_back(v.first);
            }
        } 

        if(!vertsRemaining) {
            if(DEBUG)std::cout << "No verts remaining! " << std::endl;
            return false;
        }*/

        if(DEBUG)std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        if(DEBUG)std::cout << "Verts remaining:" << std::endl;
        if(DEBUG)printVector(temp);
        begin = std::chrono::steady_clock::now();
        // Get neighbours of frontier
        //auto frontier_neighbours = diag.getConnectedSet(frontier, outputs);
        auto frontier_neighbours = get_frontier_neighbors(diag, frontier);

        if(DEBUG)std::cout << "Frontier neighbours:" << std::endl;
        if(DEBUG)printVector(frontier_neighbours);

        // Get biadjacency matrix of frontier/neighbours
        std::vector<zx::Vertex> frontier_values;
        for (const auto& [key, value] : frontier) {
            frontier_values.push_back(value);
        }

        auto adjMatrix = getAdjacencyMatrix(diag, frontier_values, frontier_neighbours);
        end = std::chrono::steady_clock::now();
        addMeasurement("extract:extractCNOT:biadjacencyMatrix", begin, end);

        if(DEBUG)std::cout << "Adjacency Matrix:" << std::endl;
        if(DEBUG)printMatrix(adjMatrix);

        // Gauss reduction on biadjacency matrix
        begin = std::chrono::steady_clock::now();
        std::vector<std::pair<zx::Qubit, zx::Qubit>> rowOperations = gaussEliminationAlt(adjMatrix);
        
        //std::set<std::pair<zx::Qubit, zx::Qubit>> unique_rowOperations(rowOperations.begin(), rowOperations.end()); // CHECK: Is this correct?
        //rowOperations.assign(unique_rowOperations.begin(), unique_rowOperations.end());

        end = std::chrono::steady_clock::now();
        addMeasurement("extract:extractCNOT:gaussElimination", begin, end);
        if(DEBUG)std::cout << "After Gauss Elimination:" << std::endl;
        if(DEBUG)printMatrix(adjMatrix);
        if(DEBUG)std::cout << "Row Operations:" << std::endl;
        if(DEBUG)printVector(rowOperations);

        std::vector<zx::Vertex> ws;
        //std::map<zx::Qubit, zx::Vertex>& ws_neighbours; // FIXME: What is this used for?
        bool singleOneRowExists = false;
        for(size_t i = 0; i < adjMatrix.size(); ++i) {
            int sum = 0;//, nonZero = 0;
            for(size_t j = 0; j < adjMatrix[i].size(); ++j) {
                sum += adjMatrix[i][j];
                /* if(adjMatrix[i][j]) {
                    nonZero = j;
                } */
            }
            if(sum == 1) {
                singleOneRowExists = true;
                break;
                // Set w to vertex corresponding to nonZero
                //int w = frontier_neighbours[nonZero];
                //ws.emplace_back(w); // FIXME: Check for duplicates (?)
            }
        }
        //std::cout << "Vector ws:" << std::endl;
        //if(DEBUG)printVector(ws);

        begin = std::chrono::steady_clock::now();
        if(!singleOneRowExists) {
            if(DEBUG)std::cout << "Ws is 0" << std::endl;
            exit(0);
            
            // Extract YZ-spiders
            for(auto v : frontier_neighbours) {
                auto v_neighbours = diag.getNeighbourVertices(v);

                //v_neighbours.erase(std::remove_if(v_neighbours.begin(), v_neighbours.end(), [diag](int x) { return diag.isInput(x); }), v_neighbours.end());

                for(auto w : v_neighbours) {
                    if(diag.isInput(w) || diag.isOutput(w)) continue; // BUG: our output...
                     auto w_neighbours = diag.getNeighbourVertices(w);

                     if(w_neighbours.size() == 1) { // BUG: Currently this can get input spiders...
                        if(DEBUG)std::cout << "Vertex with only one neighbour found: " << w << " with phase " << diag.phase(w) << std::endl;
                        if(diag.isInput(w)) {
                            if(DEBUG)std::cout << "ERROR: vertex is input!" << std::endl;
                        }

                        if(diag.phase(v).isZero()) { // Phase-gadget found 
                        //FIXME: Or PI?

                            size_t frontier_neighbour;
                            zx::Qubit q;

                            for(auto z : v_neighbours) {
                                auto it = std::find_if(frontier.begin(), frontier.end(), [z](const auto& p) { return p.second == z; });
                                if(it != frontier.end()) {
                                    frontier_neighbour = z;
                                    q = it->first;
                                    break;
                                }
                            }
                            if(frontier_neighbour) {
                                zx::pivot(diag, v, frontier_neighbour);

                                // Remove YZ-spider
                                if(DEBUG)std::cout << "Old spider " << frontier[q] << std::endl;
                                if(frontier[q] != v) {
                                    if(DEBUG)std::cout << "Old spider " << frontier[q] << " != " << v << std::endl;
                                    exit(0);
                                }
                                frontier[q] = w;

                                //frontier.erase(std::remove(frontier.begin(), frontier.end(), w), frontier.end());
                                //frontier.erase(frontier_neighbour); // FIXME: Should be qubit!

                                if(DEBUG)std::cout << "Removing YZ-spider " << v << std::endl;
                                if(DEBUG)std::cout << "New frontier spider is " << w << " on Qubit" << q << std::endl;
                            }

                            
                        }
                    }
                }
            }
            return;
        }
        end = std::chrono::steady_clock::now();
        addMeasurement("extract:extractCNOT:YZSpiders", begin, end);
        begin = std::chrono::steady_clock::now();
        // Extract CNOTs
        for(auto r : rowOperations) { // BUG: Potentially there is a bug here?
            // From row operation: r.second = r.second + r.first
            auto it = frontier.begin();
            auto control_entry = std::next(it, r.second);
            auto target_entry = std::next(it, r.first);

            auto control_qubit = control_entry->first;
            auto target_qubit = target_entry->first;
            if(DEBUG)std::cout << " Entries " << control_qubit << "|" << control_entry->second << " , " << target_qubit << "|" << target_entry->second << std::endl;

            circuit.x(target_qubit, dd::Control{(dd::Qubit) control_qubit});
            if(DEBUG)std::cout << "Added CNOT (T:" << target_qubit << ", C:" << control_qubit << ")" << std::endl;

            auto ftarg = frontier.at(control_qubit);
            auto fcont = frontier.at(target_qubit);

            // Update diag based on row operation
            for(auto v : diag.getNeighbourVertices(fcont)) {

                if( contains(outputs, v) ) {
                    continue;
                }

                if(diag.connected(ftarg, v)) {
                    diag.removeEdge(ftarg, v);
                    if(DEBUG)std::cout << "Removed edge (" << ftarg << ", " << v << ")" << std::endl;
                }
                else {
                    if(diag.isInput(v)) { // v is an input 
                        if(DEBUG)std::cout << "ASDASD " << ftarg << " vs " << r.first << std::endl;
                        auto new_v = diag.insertIdentity(fcont, target_qubit, v);
                        if(new_v) {
                            diag.addEdge(ftarg, *new_v, zx::EdgeType::Hadamard);
                            if(DEBUG)std::cout << "Added edge (" << ftarg << ", " << *new_v << ")" << std::endl;
                        }
                    }
                    else {
                        diag.addEdge(ftarg, v, zx::EdgeType::Hadamard);
                        if(DEBUG)std::cout << "Added edge (" << ftarg << ", " << v << ")" << std::endl;
                    }
                }
            }
        }
        end = std::chrono::steady_clock::now();
        addMeasurement("extract:extractCNOT:CNOTFromOperations", begin, end);

        /* for(size_t i = 0; i < ws.size(); ++i) {
            zx::Vertex v = ws_neighbours[i];
            zx::Vertex w = ws[i];
            circuit.h(diag.qubit(v)); 
            circuit.phase(diag.qubit(v), diag.phase(w).getConst().toDouble()); 
            frontier.erase(std::remove(frontier.begin(), frontier.end(), v), frontier.end());  // IMPROVE: Should really have been a map qubit->vert to begin with
            frontier.emplace_back(w);
            diag.removeVertex(v);
        } */

        /* std::vector<zx::Vertex> new_frontier;
        auto inputs = diag.getInputs();

        for(auto v : frontier) {
            auto v_neighbours = diag.getConnectedSet(frontier);
            if(v_neighbours.size() > 2 || std::find(inputs.begin(), inputs.end(), v) != inputs.end() ) {
                continue;
            }

            while(true) {

            }
        }
 */
        /* if(DEBUG)std::cout << "Frontier:" << std::endl;
        if(DEBUG)printVector(frontier);

        std::set<std::pair<Vertex, Vertex>> ws_edges = getEdgesBetweenVertices(diag, ws);
        if(DEBUG)std::cout << "ws_edges:" << std::endl;
        if(DEBUG)printVector(ws_edges);

        for(std::pair<Vertex, Vertex> e : ws_edges) {
            dd::Qubit qv = diag.qubit(e.second);
            circuit.z(diag.qubit(e.first), dd::Control{qv});
            diag.removeEdge(e.first, e.second);
        } */
        return;
    }

    

    void testExtraction(std::string filename, std::string measurementGroup) {

        if(DEBUG)std::cout << "Setting up...\n";
        qc::QuantumComputation qc{};
        std::cout << "Circuit " << filename << ":" << std::endl;
        qc.import("H:/Uni/Masterarbeit/qcec/" + filename);
        /* qc.addQubitRegister(4U);
        qc.z(0, 3_pc);
        qc.z(1, 2_pc);
        qc.h(3);
        qc.phase(1, dd::PI_2);
        qc.phase(0, dd::PI_2);
        
        qc.x(0); 
        //qc.x(1);
        qc.h(3); 
        qc.x(3, 1_pc); 
        qc.x(2, 3_pc); 
        qc.x(3, 2_pc); 
        qc.t(0);
        qc.t(1);
        qc.t(2);
        qc.tdag(3);
        qc.x(1, 0_pc);
        qc.x(3, 2_pc);
        qc.x(0, 3_pc);
        qc.x(2, 1_pc);
        qc.x(1, 0_pc);
        qc.x(3, 2_pc);
        qc.z(3, 2_pc);
        //qc.z(2, 3_pc);
        //qc.y(3, 2_pc);
        qc.tdag(0);
        qc.tdag(1);
        qc.tdag(2);
        qc.t(3);
        qc.x(1, 0_pc);
        qc.x(3, 2_pc);
        qc.s(3);
        qc.x(0, 3_pc);
        qc.h(3); */
        qc.dump("H:/Uni/Masterarbeit/qcec/original.qasm");

        if(DEBUG)std::cout << "Circuit to extract:" << std::endl;
        if(DEBUG)std::cout << qc << std::endl;
        auto begin = std::chrono::steady_clock::now();

        zx::ZXDiagram zxDiag = zx::FunctionalityConstruction::buildFunctionality(&qc);
        auto end = std::chrono::steady_clock::now();
        addMeasurement("buildFunctionality", begin, end);

        //zx::fullReduce(zxDiag);
        zxDiag.toGraphlike();
        

        std::cout << "Simplifying" << std::endl;
        begin = std::chrono::steady_clock::now();
        zx::interiorCliffordSimp(zxDiag);
        //zx::fullReduce(zxDiag);
        end = std::chrono::steady_clock::now();
        addMeasurement("interiorCliffordSimp", begin, end);

        std::cout << "Extracting" << std::endl;
        qc::QuantumComputation qc_extracted = qc::QuantumComputation(zxDiag.getNQubits());
        begin = std::chrono::steady_clock::now();
        extract(qc_extracted, zxDiag);
        end = std::chrono::steady_clock::now();
        addMeasurement("extract", begin, end);

        std::cout << "Finished Circuit" << std::endl;
        std::cout << qc_extracted << std::endl;
        //std::cout << "Circuit to extract:" << std::endl;
        //std::cout << qc << std::endl;
        qc_extracted.dump("H:/Uni/Masterarbeit/qcec/extracted.qasm");
        std::cout << "Circuit " << filename << ":" << std::endl;
        
        printMeasurements(measurementGroup, filename);

        /* ec::Configuration config{};
        config.functionality.traceThreshold    = 1e-2;
        config.execution.runAlternatingChecker = true;
        ec::EquivalenceCheckingManager ecm(qc, qc_extracted, config);
        ecm.run();
        if(DEBUG)std::cout << ecm << std::endl;
        EXPECT_EQ(ecm.equivalence(), ec::EquivalenceCriterion::Equivalent); */

    }
}



    

   



