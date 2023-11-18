#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <unordered_map>
#include "search_BWT.cpp"
#include "construct_bwt_SA.cpp"
#include "construct_output_genome.cpp"
using namespace std;

int main(){

    //create bwt and sa
    std::pair<string, vector<int>> bwt_SA = construct_bwt_SA();
    string bwt = bwt_SA.first;
    vector<int> SA = bwt_SA.second;

    //get the genome and reads
    std::string genome;
    vector<string> reads;
    try {
        genome = read_genome("genome.txt");
        reads = read_reads("reads.txt");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    //indexing
    // Start indexing timer
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<int>, std::vector<std::vector<int>>> c_occ = index_bwt(bwt);
    // Stop indexing timer
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "indexing time: " << duration.count()/1000.0 << " seconds" << std::endl;

    vector<int> c = c_occ.first;
    vector<vector<int>> occ = c_occ.second;

    vector<pair<int, string>> alignments;
    
    //aligning
    // Start aligning timer
    start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<reads.size(); i++)
    {
        string read = reads[i];
        // cout<<"Processing read "<<i<<" out of "<<reads.size()<<'\n';
        vector<int> hits = search(bwt, read, SA, c, occ);
        for(auto hit : hits){
            pair<int, string> alignment = {hit, read};
            alignments.push_back(alignment);
        }
    }
    // Stop aligning timer
    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Aligning time: " << duration.count()/1000.0 << " seconds" << std::endl;

    writeToFile("assembled_genome.txt", alignments);
    return 0;
}