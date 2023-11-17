#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
#include "libsais.h"

std::string read_genome(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }
    std::string genome((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    return genome;
}

std::vector<std::string> read_reads(const std::string& path) {
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

int construct_bwt_SA() {
    std::string file_path = "./genome.txt";

    std::ifstream file(file_path);
    if (!file) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return 1;
    }

    std::string nucleotides;
    std::string line;
    while (file >> line) {
        nucleotides += line;
    }
    file.close();

    std::vector<uint8_t> uint8_nucleotides(nucleotides.begin(), nucleotides.end());
    std::vector<int32_t> suffix_array(uint8_nucleotides.size(), 0);
    std::vector<uint8_t> bwt_output(uint8_nucleotides.size(), 0);

    auto start_sa = std::chrono::high_resolution_clock::now();

    int32_t status_sa = libsais(
        uint8_nucleotides.data(),
        suffix_array.data(),
        static_cast<int32_t>(uint8_nucleotides.size()),
        0,
        NULL
    );

    auto end_sa = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_sa = end_sa - start_sa;

    if (status_sa != 0) {
        std::cerr << "Error constructing suffix array" << std::endl;
        return 1;
    }

    // Saving the suffix array to SA_Array.txt
    std::ofstream output_sa("./SA.txt");
    output_sa<<nucleotides.size()<<std::endl;
    for (size_t i = 0; i < suffix_array.size(); ++i) {
        output_sa << suffix_array[i] << std::endl;
    }
    output_sa.close(); // Close the file

    auto start_bwt = std::chrono::high_resolution_clock::now();

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
        return 1;
    }

    std::ofstream output_bwt("./bwt.txt");
    for (size_t i = 0; i < bwt_output.size(); ++i) {
        if (i == static_cast<size_t>(primary_index)) {
            output_bwt << '$';
        }
        output_bwt << static_cast<char>(bwt_output[i]);
    }
    output_bwt.close(); // Close the file

    std::cout << "Time taken for Suffix Array construction: " << duration_sa.count() << " milliseconds" << '\n';
    std::cout << "Time taken for BWT construction: " << duration_bwt.count() << " milliseconds" << '\n';

    return 0;
}
