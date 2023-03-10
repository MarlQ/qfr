#include "zx/Extraction.hpp"
#include "zx/Measurement.hpp"

#include "Definitions.hpp"
#include "Rational.hpp"
#include "Utils.hpp"
#include "ZXDiagram.hpp"
#include "Simplify.hpp"
#include "dd/Control.hpp"
#include "zx/FunctionalityConstruction.hpp"
#include "CircuitOptimizer.hpp"
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
#include <omp.h>

namespace zx {

    class Extractor {
        public:
            Extractor(qc::QuantumComputation& circuit, ZXDiagram& diag, bool parallelize=false) 
            : circuit(circuit), diag(diag), parallelize(parallelize), inputs(diag.getInputs()), outputs(diag.getOutputs())
            {
                frontier = initFrontier();
            };

            Extractor* other_extractor;

            void extract() {
                // TODO: ...



                // Check whether frontier vertices are adjacent between the two extractors
                // Extractor 2 removes the overlapping adjcacent vertices from its scope
                if(parallelize) {
                    #pragma omp barrier

                    #pragma omp critical
                    {
                        // Determine overlapping elements between the two Extractor objects
                        std::vector<size_t> intersection;
                        std::set_intersection(frontier_neighbours.begin(), frontier_neighbours.end(),
                                            other_extractor->frontier_neighbours.begin(), other_extractor->frontier_neighbours.end(),
                                            std::back_inserter(intersection));

                        // Remove overlapping elements from one of the Extractor objects
                        if (intersection.size() > 0) {
                            // Remove overlapping elements from one of the Extractor objects
                            if (omp_get_thread_num() == 0) {
                                /* frontier_neighbours.erase(
                                    std::remove_if(extractor1.frontier_neighbours.begin(), extractor1.frontier_neighbours.end(),
                                                [&intersection](size_t i) { return std::find(intersection.begin(), intersection.end(), i) != intersection.end(); }),
                                    extractor1.frontier_neighbours.end()); */
                            }
                        }
                    }

                    #pragma omp barrier
                }
                
                // TODO: ...


            }

        private:
            std::vector<size_t> frontier_neighbours;
            qc::QuantumComputation& circuit;
            ZXDiagram& diag;
            std::vector<size_t> inputs;
            std::vector<size_t> outputs;
            std::map<zx::Qubit, zx::Vertex> frontier;
            bool parallelize;
            

            std::map<zx::Qubit, zx::Vertex> initFrontier() {
                std::map<zx::Qubit, zx::Vertex> frontier;
                for(size_t i = 0; i < outputs.size(); ++i) {
                    auto v = diag.getNeighbourVertices(outputs[i])[0];
                    if( !diag.isInput(v) ) {
                        frontier[i] = v;
                    }
                }

                return frontier;
            }
    };

    #define DEBUG false

    void testExtraction(std::string circuitName, std::string measurementGroup) {
        bool parallelize = true;
        Measurement measurement;

        if(DEBUG)std::cout << "Setting up...\n";
        qc::QuantumComputation qc{};
        std::cout << "Circuit " << circuitName << ":" << std::endl;
        qc.import("H:/Uni/Masterarbeit/qcec/" + circuitName);

        qc.dump("H:/Uni/Masterarbeit/qcec/original.qasm");

        if(DEBUG)std::cout << "Circuit to extract:" << std::endl;
        if(DEBUG)std::cout << qc << std::endl;
        auto begin = std::chrono::steady_clock::now();

        zx::ZXDiagram zxDiag = zx::FunctionalityConstruction::buildFunctionality(&qc);
        auto end = std::chrono::steady_clock::now();
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

        Extractor extractor1(qc_extracted, zxDiag);


        if(parallelize) {
            extractor1.extract();
        }
        else {
            zx::ZXDiagram zxDiag_reversed = zxDiag.adjoint();
            qc::QuantumComputation qc_extracted_2 = qc::QuantumComputation(zxDiag_reversed.getNQubits());

            Extractor extractor2(qc_extracted_2, zxDiag_reversed);

            extractor1.other_extractor = &extractor2;
            extractor2.other_extractor = &extractor1;

            #pragma omp parallel num_threads(2)
            {
                if(omp_get_thread_num() == 0) {
                    extractor1.extract();
                }
                else {
                    extractor2.extract();
                }
            }

        }

        





        end = std::chrono::steady_clock::now();
        measurement.addMeasurement("extract", begin, end);

        std::cout << "Finished Circuit" << std::endl;
        std::cout << qc_extracted << std::endl;
        //std::cout << "Circuit to extract:" << std::endl;
        //std::cout << qc << std::endl;
        qc_extracted.dump("H:/Uni/Masterarbeit/pyzx/thesis/extracted.qasm");
        std::cout << "Circuit " << circuitName << ":" << std::endl;
        
        measurement.printMeasurements(measurementGroup, circuitName);

        

    }

}