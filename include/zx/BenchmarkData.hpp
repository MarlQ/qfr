#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

class BenchmarkData {
public:
    std::string measurement_group = "MEASUREMENT";
    std::string circuit_name = "CIRCUIT";
    double time_total = 0;
    double time_parallel = 0;
    int parallel_iterations = 0;
    int total_iterations = 0;
    int iteration_diff = 0;

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

    int num_extr_par_cnot = 0;
    int num_extr_par_cz = 0;
    int num_extr_par_fp = 0;

    int num_extr_seq_cnot = 0;
    int num_extr_seq_cz = 0;
    int num_extr_seq_fp = 0;

    int failedCnots = 0;

    int num_gates_cnot = 0;
    int num_gates_cz = 0;
    int num_gates_phase = 0;
    int num_gates_h = 0;
    int num_gates_swap = 0;

        // Merge data from another BenchmarkData object and average the fields
    void averageData(const BenchmarkData& other) {
        time_total = (time_total + other.time_total);
        time_parallel = time_parallel + other.time_parallel;
        parallel_iterations = (parallel_iterations + other.parallel_iterations);
        total_iterations = (total_iterations + other.total_iterations);
        iteration_diff = iteration_diff + other.iteration_diff;
        time_extr_par_cnot = std::max(time_extr_par_cnot, other.time_extr_par_cnot);
        time_extr_par_cz = std::max(time_extr_par_cz, other.time_extr_par_cz);
        time_extr_par_fp = std::max(time_extr_par_fp, other.time_extr_par_fp);
        time_extr_seq_cnot = (time_extr_seq_cnot+ other.time_extr_seq_cnot);
        time_extr_seq_cz = (time_extr_seq_cz+ other.time_extr_seq_cz);
        time_extr_seq_fp = (time_extr_seq_fp+ other.time_extr_seq_fp);
        time_cnot_failed_extraction = std::max(time_cnot_failed_extraction,other.time_cnot_failed_extraction);
        time_cnot_gauss = (time_cnot_gauss+ other.time_cnot_gauss);
        time_cnot_biadj = (time_cnot_biadj+ other.time_cnot_biadj);
        time_cnot_optimal = (time_cnot_optimal+ other.time_cnot_optimal);
        time_cnot_neighbors = (time_cnot_neighbors+ other.time_cnot_neighbors);
        num_extr_par_cnot = (num_extr_par_cnot + other.num_extr_par_cnot);
        num_extr_par_cz = (num_extr_par_cz + other.num_extr_par_cz);
        num_extr_par_fp = (num_extr_par_fp + other.num_extr_par_fp);
        num_extr_seq_cnot = (num_extr_seq_cnot + other.num_extr_seq_cnot);
        num_extr_seq_cz = (num_extr_seq_cz + other.num_extr_seq_cz);
        num_extr_seq_fp = (num_extr_seq_fp + other.num_extr_seq_fp);
        failedCnots = (failedCnots + other.failedCnots);
        num_gates_cnot = (num_gates_cnot + other.num_gates_cnot);
        num_gates_cz = (num_gates_cz + other.num_gates_cz);
        num_gates_phase = (num_gates_phase + other.num_gates_phase);
        num_gates_h = (num_gates_h + other.num_gates_h);
        num_gates_swap = (num_gates_swap + other.num_gates_swap);
    }

    void mergeData(const BenchmarkData& other) {
        dataAggregates++;
        //time_total = (time_total+ other.time_total);
        parallel_iterations = (parallel_iterations + other.parallel_iterations);
        total_iterations = (total_iterations + other.total_iterations);
        iteration_diff = iteration_diff + other.iteration_diff;
        time_extr_par_cnot = (time_extr_par_cnot+ other.time_extr_par_cnot);
        time_extr_par_cz = (time_extr_par_cz+ other.time_extr_par_cz);
        time_extr_par_fp = (time_extr_par_fp+ other.time_extr_par_fp);
        time_extr_seq_cnot = (time_extr_seq_cnot+ other.time_extr_seq_cnot);
        time_extr_seq_cz = (time_extr_seq_cz+ other.time_extr_seq_cz);
        time_extr_seq_fp = (time_extr_seq_fp+ other.time_extr_seq_fp);
        time_cnot_failed_extraction = (time_cnot_failed_extraction+ other.time_cnot_failed_extraction);
        time_cnot_gauss = (time_cnot_gauss+ other.time_cnot_gauss);
        time_cnot_biadj = (time_cnot_biadj+ other.time_cnot_biadj);
        time_cnot_optimal = (time_cnot_optimal+ other.time_cnot_optimal);
        time_cnot_neighbors = (time_cnot_neighbors+ other.time_cnot_neighbors);
        num_extr_par_cnot = (num_extr_par_cnot + other.num_extr_par_cnot);
        num_extr_par_cz = (num_extr_par_cz + other.num_extr_par_cz);
        num_extr_par_fp = (num_extr_par_fp + other.num_extr_par_fp);
        num_extr_seq_cnot = (num_extr_seq_cnot + other.num_extr_seq_cnot);
        num_extr_seq_cz = (num_extr_seq_cz + other.num_extr_seq_cz);
        num_extr_seq_fp = (num_extr_seq_fp + other.num_extr_seq_fp);
        failedCnots = (failedCnots + other.failedCnots);
        num_gates_cnot = (num_gates_cnot + other.num_gates_cnot);
        num_gates_cz = (num_gates_cz + other.num_gates_cz);
        num_gates_phase = (num_gates_phase + other.num_gates_phase);
        num_gates_h = (num_gates_h + other.num_gates_h);
        num_gates_swap = (num_gates_swap + other.num_gates_swap);
        time_totals.emplace_back(other.time_total);
        time_parallel_totals.emplace_back(other.time_parallel);

        //time_total = std::min(time_total, other.time_total);
/*         parallel_iterations = std::min(parallel_iterations , other.parallel_iterations);
        total_iterations = std::min(total_iterations , other.total_iterations);
        time_extr_par_cnot = std::min(time_extr_par_cnot, other.time_extr_par_cnot);
        time_extr_par_cz = std::min(time_extr_par_cz, other.time_extr_par_cz);
        time_extr_par_fp = std::min(time_extr_par_fp, other.time_extr_par_fp);
        time_extr_seq_cnot = std::min(time_extr_seq_cnot, other.time_extr_seq_cnot);
        time_extr_seq_cz = std::min(time_extr_seq_cz, other.time_extr_seq_cz);
        time_extr_seq_fp = std::min(time_extr_seq_fp, other.time_extr_seq_fp);
        time_cnot_failed_extraction = std::min(time_cnot_failed_extraction, other.time_cnot_failed_extraction);
        time_cnot_gauss = std::min(time_cnot_gauss, other.time_cnot_gauss);
        time_cnot_biadj = std::min(time_cnot_biadj, other.time_cnot_biadj);
        time_cnot_optimal = std::min(time_cnot_optimal, other.time_cnot_optimal);
        time_cnot_neighbors = std::min(time_cnot_neighbors, other.time_cnot_neighbors);
        num_extr_par_cnot = std::min(num_extr_par_cnot , other.num_extr_par_cnot);
        num_extr_par_cz = std::min(num_extr_par_cz , other.num_extr_par_cz);
        num_extr_par_fp = std::min(num_extr_par_fp , other.num_extr_par_fp);
        num_extr_seq_cnot = std::min(num_extr_seq_cnot , other.num_extr_seq_cnot);
        num_extr_seq_cz = std::min(num_extr_seq_cz , other.num_extr_seq_cz);
        num_extr_seq_fp = std::min(num_extr_seq_fp , other.num_extr_seq_fp);
        failedCnots = std::min(failedCnots , other.failedCnots);
        num_gates_cnot = std::min(num_gates_cnot , other.num_gates_cnot);
        num_gates_cz = std::min(num_gates_cz , other.num_gates_cz);
        num_gates_phase = std::min(num_gates_phase , other.num_gates_phase);
        num_gates_h = std::min(num_gates_h , other.num_gates_h);
        num_gates_swap = std::min(num_gates_swap , other.num_gates_swap); */
        
    }

    void printStatistics() { 
        // Print the statistics
/*         std::cout << "-----------------------------------------------" << std::endl;

        std::cout << "// Time for extraction operations" << std::endl;
        std::cout << "time_extr_par_cnot = " << time_extr_par_cnot  << std::endl;
        std::cout << "time_extr_par_cz = " << time_extr_par_cz  << std::endl;
        std::cout << "time_extr_par_fp = " << time_extr_par_fp  << std::endl;
        std::cout << std::endl;

        std::cout << "time_extr_seq_cnot = " << time_extr_seq_cnot  << std::endl;
        std::cout << "time_extr_seq_cz = " << time_extr_seq_cz  << std::endl;
        std::cout << "time_extr_seq_fp = " << time_extr_seq_fp  << std::endl;
        std::cout << std::endl;

        std::cout << "time_cnot_failed_extraction = " << time_cnot_failed_extraction  << std::endl;
        std::cout << std::endl;

        std::cout << "// Number of extraction operations" << std::endl;
        std::cout << "num_extr_par_cnot = " << num_extr_par_cnot  << std::endl;
        std::cout << "num_extr_par_cz = " << num_extr_par_cz  << std::endl;
        std::cout << "num_extr_par_fp = " << num_extr_par_fp  << std::endl;
        std::cout << std::endl;

        std::cout << "num_extr_seq_cnot = " << num_extr_seq_cnot  << std::endl;
        std::cout << "num_extr_seq_cz = " << num_extr_seq_cz  << std::endl;
        std::cout << "num_extr_seq_fp = " << num_extr_seq_fp  << std::endl;
        std::cout << std::endl;
        
        std::cout << "failedCnots = " << failedCnots  << std::endl;
        std::cout << std::endl;
        

        std::cout << "// Number of gates created during extraction" << std::endl;
        std::cout << "num_gates_cnot = " << num_gates_cnot  << std::endl;
        std::cout << "num_gates_cz = " << num_gates_cz  << std::endl;
        std::cout << "num_gates_phase = " << num_gates_phase  << std::endl;
        std::cout << "num_gates_h = " << num_gates_h  << std::endl;
        std::cout << "num_gates_swap = " << num_gates_swap  << std::endl; */


        // Write the output to a CSV file
        std::ostringstream output;

        output << measurement_group;
        output << "," << circuit_name;
        output << "," << time_total;
        output << "," << time_parallel;
        output << "," << parallel_iterations;
        output << "," << total_iterations;
        output << "," << (total_iterations > 0 ?  ((double)parallel_iterations / (double)total_iterations) : 0);
        output << "," << iteration_diff;
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

    void finish() {
        time_total = median(time_totals);
        time_parallel = median(time_parallel_totals);
        parallel_iterations /=dataAggregates;
        total_iterations /=dataAggregates;
        iteration_diff /= dataAggregates;

        time_extr_par_cnot /=dataAggregates;
        time_extr_par_cz /=dataAggregates;
        time_extr_par_fp /=dataAggregates;

        time_extr_seq_cnot /=dataAggregates;
        time_extr_seq_cz /=dataAggregates;
        time_extr_seq_fp /=dataAggregates;

        time_cnot_failed_extraction /=dataAggregates;

        time_cnot_gauss /=dataAggregates;
        time_cnot_biadj /=dataAggregates;
        time_cnot_optimal /=dataAggregates;
        time_cnot_neighbors /=dataAggregates;

        num_extr_par_cnot /=dataAggregates;
        num_extr_par_cz /=dataAggregates;
        num_extr_par_fp /=dataAggregates;

        num_extr_seq_cnot /=dataAggregates;
        num_extr_seq_cz /=dataAggregates;
        num_extr_seq_fp /=dataAggregates;

        failedCnots /=dataAggregates;

        num_gates_cnot /=dataAggregates;
        num_gates_cz /=dataAggregates;
        num_gates_phase /=dataAggregates;
        num_gates_h /=dataAggregates;
        num_gates_swap /=dataAggregates;

        printStatistics();
        // Reset
        time_total = 0;
        parallel_iterations = 0;
        total_iterations = 0;
        iteration_diff = 0;

        time_extr_par_cnot = 0;
        time_extr_par_cz = 0;
        time_extr_par_fp = 0;

        time_extr_seq_cnot = 0;
        time_extr_seq_cz = 0;
        time_extr_seq_fp = 0;

        time_cnot_failed_extraction = 0;

        time_cnot_gauss = 0;
        time_cnot_biadj = 0;
        time_cnot_optimal = 0;
        time_cnot_neighbors = 0;

        num_extr_par_cnot = 0;
        num_extr_par_cz = 0;
        num_extr_par_fp = 0;

        num_extr_seq_cnot = 0;
        num_extr_seq_cz = 0;
        num_extr_seq_fp = 0;

        failedCnots = 0;

        num_gates_cnot = 0;
        num_gates_cz = 0;
        num_gates_phase = 0;
        num_gates_h = 0;
        num_gates_swap = 0;

        dataAggregates = 0;

        time_totals.clear();
        time_parallel_totals.clear();
    }

private:
    int dataAggregates = 0;
    double average(double a, double b) {
        return (a+b)/2.0;
    }

    int average(int a, int b) {
        return (a+b)/2;
    }

/*     int median(std::vector<int> &v) {
        size_t n = v.size() / 2;
        std::nth_element(v.begin(), v.begin()+n, v.end());
        return v[n];
    } */

    double median(std::vector<double> &v) {
        size_t n = v.size() / 2;
        std::nth_element(v.begin(), v.begin()+n, v.end());
        return v[n];
    }

    std::vector<double> time_totals;
    std::vector<double> time_parallel_totals;
};
