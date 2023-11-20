#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
#include "libsais.h"
#include <algorithm> 
#include <cctype>

using namespace std;

//read genome from file
string read_genome(const string& path) {
    ifstream file(path);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + path);
    }
    // Read entire file content into a string
    string genome((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
    // Convert to uppercase
    std::transform(genome.begin(), genome.end(), genome.begin(), ::toupper);

    return genome;
}

//function to construct bwt and SA and return them
pair<string, vector<int>> construct_bwt_SA () {
    string file_path = "./reference.txt";

    ifstream file(file_path);
    if (!file) {
        cerr << "Error opening file: " << file_path << endl;
        return {};
    }

    string nucleotides;
    string line;
    while (file >> line) {
        nucleotides += line;
    }
    
    std::transform(nucleotides.begin(), nucleotides.end(), nucleotides.begin(), ::toupper);
    file.close();

    vector<uint8_t> uint8_nucleotides(nucleotides.begin(), nucleotides.end());
    vector<int> suffix_array(uint8_nucleotides.size(), 0);
    vector<uint8_t> bwt_output(uint8_nucleotides.size(), 0);

    auto start_sa = chrono::high_resolution_clock::now();
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
    auto end_sa = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_sa = end_sa - start_sa;

    if (status_sa != 0) {
        cerr << "Error constructing suffix array" << endl;
        return {};
    }

    auto start_bwt = chrono::high_resolution_clock::now();
    //bwt construction
    int32_t primary_index = libsais_bwt(
        uint8_nucleotides.data(),
        bwt_output.data(),
        suffix_array.data(),
        static_cast<int32_t>(uint8_nucleotides.size()),
        0,
        NULL
    );

    auto end_bwt = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_bwt = end_bwt - start_bwt;

    if (primary_index < 0) {
        cerr << "Error constructing BWT" << endl;
        return {};
    }

    cout << "Time taken for Suffix Array construction: " << duration_sa.count() << " milliseconds" << '\n';
    cout << "Time taken for BWT construction: " << duration_bwt.count() << " milliseconds" << '\n';

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
