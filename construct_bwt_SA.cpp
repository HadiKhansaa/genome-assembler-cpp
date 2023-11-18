#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
#include "libsais.h"
using namespace std;

std::string read_genome(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }
    std::string genome((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    return genome;
}

std::vector<std::string> read_reads(const std::string& path) {
    //start timer
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }

    std::vector<std::string> reads;
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            reads.push_back(line);
        }
    }
    //end timer
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "Time taken to read reads: " << duration.count()/1000.0 << " seconds" << std::endl;
    return reads;
}

std::string read_bwt(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }
    std::string bwt((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    return bwt;
}

std::vector<int> read_sa(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }

    std::vector<int> suffixArray;
    int value;
    while (file >> value) {
        suffixArray.push_back(value);
    }

    return suffixArray;
}

//function to construct bwt and SA and return them
pair<string, vector<int>> construct_bwt_SA () {
    std::string file_path = "./genome.txt";

    std::ifstream file(file_path);
    if (!file) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return {};
    }

    std::string nucleotides;
    std::string line;
    while (file >> line) {
        nucleotides += line;
    }
    file.close();

    std::vector<uint8_t> uint8_nucleotides(nucleotides.begin(), nucleotides.end());
    std::vector<int> suffix_array(uint8_nucleotides.size(), 0);
    std::vector<uint8_t> bwt_output(uint8_nucleotides.size(), 0);

    auto start_sa = std::chrono::high_resolution_clock::now();
    //suffix array construction
    int32_t status_sa = libsais(
        uint8_nucleotides.data(),
        suffix_array.data(),
        static_cast<int32_t>(uint8_nucleotides.size()),
        0,
        NULL
    );

    vector<int> SA;
    SA.push_back(nucleotides.size());
    for (size_t i = 0; i < suffix_array.size(); ++i) {
        SA.push_back(suffix_array[i]);
    }
    auto end_sa = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_sa = end_sa - start_sa;

    if (status_sa != 0) {
        std::cerr << "Error constructing suffix array" << std::endl;
        return {};
    }

    auto start_bwt = std::chrono::high_resolution_clock::now();
    //bwt construction
    int32_t primary_index = libsais_bwt(
        uint8_nucleotides.data(),
        bwt_output.data(),
        suffix_array.data(),
        static_cast<int32_t>(uint8_nucleotides.size()),
        0,
        NULL
    );

    auto end_bwt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_bwt = end_bwt - start_bwt;

    if (primary_index < 0) {
        std::cerr << "Error constructing BWT" << std::endl;
        return {};
    }

    std::cout << "Time taken for Suffix Array construction: " << duration_sa.count() << " milliseconds" << '\n';
    std::cout << "Time taken for BWT construction: " << duration_bwt.count() << " milliseconds" << '\n';

    string bwt_string = "";
    for (size_t i = 0; i < bwt_output.size(); ++i) {
        if (i == static_cast<size_t>(primary_index)) {
            bwt_string += '$';
        }
        bwt_string += static_cast<char>(bwt_output[i]);
    }
    
    //return bwt and SA
    return {bwt_string, SA};
}
