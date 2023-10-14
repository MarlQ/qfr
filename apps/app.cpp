/*
 * This file is part of MQT QCEC library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de//research/quantum/ for more information.
 */

#include "QuantumComputation.hpp"
#include "zx/ExtractorParallel.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <locale>
#include <random>
#include <set>
#include <string>

void show_usage(const std::string& name) {
    std::cerr << "Usage: " << name << "<PATH_INPUT_FILE> <PATH_TO_OUTPUT_FILE> (--remove_gates X)" << std::endl;
    std::cerr << "Supported input file formats:" << std::endl;
    std::cerr << "  .real                       " << std::endl;
    std::cerr << "  .qasm                       " << std::endl;
    std::cerr << "Supported output file formats:" << std::endl;
    std::cerr << "  .qasm                       " << std::endl;
    std::cerr << "  .py (qiskit)                " << std::endl;
    std::cerr << "If '--remove_gates X' is specified, X gates are randomly removed" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "Starting App" << std::endl;

    std::vector<std::string> circuits = {
        "circuits\\large\\adder_n28\\adder_n28_transpiled.qasm",
        "circuits\\large\\adder_n64\\adder_n64_transpiled.qasm",
        "circuits\\large\\adder_n118\\adder_n118_transpiled.qasm",
        "circuits\\pyzx\\barenco_tof_5.qasm",
        "circuits\\pyzx\\barenco_tof_10.qasm",
        "circuits\\large\\bv_n30\\bv_n30_transpiled.qasm",
        "circuits\\large\\bv_n70\\bv_n70_transpiled.qasm",
        "circuits\\large\\cat_n35\\cat_n35_transpiled.qasm",
        "circuits\\large\\cat_n65\\cat_n65_transpiled.qasm",
        "circuits\\pyzx\\gf2^10_mult.qasm",
        "circuits\\pyzx\\gf2^16_mult.qasm",
        "circuits\\pyzx\\gf2^32_mult.qasm", 
        "circuits\\large\\ghz_n40\\ghz_n40_transpiled.qasm",
        "circuits\\large\\ghz_n78\\ghz_n78_transpiled.qasm",
        "circuits\\large\\ghz_n127\\ghz_n127_transpiled.qasm",
        "circuits\\pyzx\\ham15-low.qasm",
        "circuits\\pyzx\\ham15-med.qasm",
        "circuits\\pyzx\\ham15-high.qasm",
        "circuits\\small\\hhl_n7\\hhl_n7_transpiled.qasm",
        "circuits\\small\\hhl_n10\\hhl_n10_transpiled.qasm",
        "circuits\\pyzx\\hwb6.qasm",
        "circuits\\pyzx\\hwb8.qasm",
        "circuits\\pyzx\\hwb10.qasm",
        "circuits\\pyzx\\hwb12.qasm",
        "circuits\\pyzx\\mod_adder_1024.qasm",
        "circuits\\large\\multiplier_n45\\multiplier_n45_transpiled.qasm", 
        "circuits\\large\\multiplier_n75\\multiplier_n75_transpiled.qasm",
        "circuits\\large\\qft_n29\\qft_n29_transpiled.qasm",
        "circuits\\large\\qft_n63\\qft_n63_transpiled.qasm",
        "circuits\\medium\\qram_n20\\qram_n20_transpiled.qasm",
        "circuits\\large\\qugan_n39\\qugan_n39_transpiled.qasm",
        "circuits\\large\\qugan_n71\\qugan_n71_transpiled.qasm",
        "circuits\\large\\qugan_n111\\qugan_n111_transpiled.qasm",
        "circuits\\pyzx\\vbe_adder_3.qasm",
        "circuits\\small\\vqe_uccsd_n4\\vqe_uccsd_n4_transpiled.qasm",
        "circuits\\small\\vqe_uccsd_n6\\vqe_uccsd_n6_transpiled.qasm",
        "circuits\\small\\vqe_uccsd_n8\\vqe_uccsd_n8_transpiled.qasm",
        "circuits\\large\\wstate_n36\\wstate_n36_transpiled.qasm",
        "circuits\\large\\wstate_n76\\wstate_n76_transpiled.qasm",
        "circuits\\large\\wstate_n118\\wstate_n118_transpiled.qasm",
    };
    bool benchmark = true;
    int benchmarkIterations = 1;
    std::string benchmarkName = "B6";
    if(benchmark) {
        zx::ExtractorConfig config;
        config.perm_optimization = true;
        int counter = 1;
        //double average = 0;
        for(std::string circuit : circuits) {
            std::cout << "Benchmark " << counter << " / " << circuits.size() << " (" << circuit <<")" << std::endl;
/* 
             std::unique_ptr<qc::QuantumComputation> qc = std::make_unique<qc::QuantumComputation>();
            qc->import("H:/Uni/Masterarbeit/qcec/" + circuit);
            zx::ZXDiagram zxDiag = zx::FunctionalityConstruction::buildFunctionality(qc.get());
            zxDiag.toGraphlike();
            int vertCountBefore = zxDiag.getNVertices();
            zx::interiorCliffordSimp(zxDiag);
            int vertCountAfter = zxDiag.getNVertices();
            double percentage = ( ((double) vertCountBefore) -  ((double) vertCountAfter) ) / ((double) vertCountBefore) * 100;
            average += percentage;
            std::cout << "Vertices: " << vertCountAfter << " / " << vertCountBefore << "(" << percentage << ")" << std::endl; 
            continue;*/
/*             std::string filename = "H:\\Uni\\Masterarbeit\\pyzx\\thesis\\" + circuit + ".json";
            std::cout << "Writing to " << filename << std::endl;

            zxDiag.toJSON(filename); */
            BenchmarkData averageData;
            averageData.circuit_name = circuit;

            std::cout << "Testing sequential version" << std::endl;
            averageData.measurement_group = benchmarkName + "_seq";
            for(int i = 0; i < benchmarkIterations; i++) {
                BenchmarkData executionData = zx::testParallelExtraction(circuit, benchmarkName + "_seq", false, config);
                averageData.mergeData(executionData);
            }
            averageData.finish();
            
            std::cout << "Testing parallel version" << std::endl;
            config.parallel_frontier_processing = false;
            averageData.measurement_group = benchmarkName + "_par";
            for(int i = 0; i < benchmarkIterations; i++) { 
                BenchmarkData executionData = zx::testParallelExtraction(circuit, benchmarkName + "_par", true, config);
                averageData.mergeData(executionData);
            }
            averageData.finish();
 
/*             std::cout << "Testing parallel version with claimed FP" << std::endl;
            config.parallel_frontier_processing = true;
            averageData.measurement_group = benchmarkName + "_parFP";
            for(int i = 0; i < benchmarkIterations; i++) { 
                std::cout << "I " << i << std::endl;
                BenchmarkData executionData = zx::testParallelExtraction(circuit, benchmarkName + "_parFP", true, config);
                averageData.mergeData(executionData);
            }
            averageData.finish(); */
            counter++;
        }
        return 0;
    }

/*     if(argc > 3) {
        bool parallelization = strcmp(argv[3], "true") == 0 || strcmp(argv[3], "1") == 0;
        std::cout << "Parallelization: " << parallelization << " | " << argv[3] << std::endl;

        if(argc > 4) {
            zx::testParallelExtraction(argv[1], argv[2], parallelization, true, std::stoi(argv[4]));
        }
        else {
            zx::testParallelExtraction(argv[1], argv[2], parallelization);   
        }
    }
    else zx::testParallelExtraction(); */
    return 0;


    if (argc != 3 && argc != 5) {
        show_usage(argv[0]);
        return 1;
    }

    // get filenames
    std::string infile  = argv[1];
    std::string outfile = argv[2];

    // get file format
    qc::Format  informat;
    size_t      dot       = infile.find_last_of('.');
    std::string extension = infile.substr(dot + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](const unsigned char c) { return ::tolower(c); });
    if (extension == "real") {
        informat = qc::Real;
    } else if (extension == "qasm") {
        informat = qc::OpenQASM;
    } else {
        show_usage(argv[0]);
        return 1;
    }

    qc::Format outformat;
    dot       = outfile.find_last_of('.');
    extension = outfile.substr(dot + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](const unsigned char c) { return ::tolower(c); });
    if (extension == "py") {
        outformat = qc::Qiskit;
    } else if (extension == "qasm") {
        outformat = qc::OpenQASM;
    } else {
        show_usage(argv[0]);
        return 1;
    }

    // read circuit
    qc::QuantumComputation qc;
    qc.import(infile, informat);

    if (argc > 3) {
        unsigned long long gates_to_remove = std::stoull(argv[4]);

        std::array<std::mt19937_64::result_type, std::mt19937_64::state_size> random_data{};
        std::random_device                                                    rd;
        std::generate(begin(random_data), end(random_data), [&]() { return rd(); });
        std::seed_seq                                     seeds(begin(random_data), end(random_data));
        std::mt19937_64                                   mt(seeds);
        std::uniform_int_distribution<unsigned long long> distribution(0, qc.getNops() - 1);
        std::function<unsigned long long()>               rng = [&]() { return distribution(mt); };

        std::set<unsigned long long> already_removed{};

        for (unsigned long long j = 0; j < gates_to_remove; ++j) {
            auto gate_to_remove = rng() % qc.getNops();
            while (already_removed.count(gate_to_remove)) {
                gate_to_remove = rng() % qc.getNops();
            }
            already_removed.insert(gate_to_remove);
            auto it = qc.begin();
            if (it == qc.end()) continue;
            std::advance(it, gate_to_remove);

            qc.erase(it);
        }
    }

    qc.dump(outfile, outformat);

    return 0;
}
